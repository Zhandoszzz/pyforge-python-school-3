from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.providers.postgres.hooks.postgres import PostgresHook
from airflow.providers.amazon.aws.hooks.s3 import S3Hook
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os
from datetime import datetime


def extract_data():
    pg_hook = PostgresHook(postgres_conn_id='postgres_default')
    connection = pg_hook.get_conn()
    cursor = connection.cursor()

    today = datetime.today().strftime('%Y-%m-%d')

    sql_query = f"""
    SELECT smiles, name
    FROM drugs
    WHERE date_column = '{today}';
    """

    cursor.execute(sql_query)
    data = cursor.fetchall()

    df = pd.DataFrame(data, columns=['smiles', 'name'])

    return df.to_dict()


def transform_data(**kwargs):
    ti = kwargs['ti']
    data_dict = ti.xcom_pull(task_ids='extract_data')

    df = pd.DataFrame.from_dict(data_dict)

    mols = [Chem.MolFromSmiles(smiles) for smiles in df['smiles']]
    df['MolecularWeight'] = [Descriptors.MolWt(mol) for mol in mols]
    df['logP'] = [Descriptors.MolLogP(mol) for mol in mols]
    df['TPSA'] = [Descriptors.TPSA(mol) for mol in mols]
    df['H_Donors'] = [Descriptors.NumHDonors(mol) for mol in mols]
    df['H_Acceptors'] = [Descriptors.NumHAcceptors(mol) for mol in mols]

    df['Lipinski_pass'] = (
        (df['MolecularWeight'] < 500) &
        (df['logP'] < 5) &
        (df['TPSA'] < 140) &
        (df['H_Donors'] <= 5) &
        (df['H_Acceptors'] <= 10)
    )

    return df.to_dict()


def save_to_s3(**kwargs):
    ti = kwargs['ti']
    data_dict = ti.xcom_pull(task_ids='transform_data')

    df = pd.DataFrame.from_dict(data_dict)

    file_name = f"transformed_data_{datetime.today().strftime('%Y-%m-%d')}.xlsx"
    df.to_excel(file_name, index=False)

    s3_hook = S3Hook(aws_conn_id='aws_default')
    s3_bucket = 'hw-bucket-yep1'
    s3_key = f"{file_name}"

    s3_hook.load_file(filename=file_name, key=s3_key, bucket_name=s3_bucket, replace=True)
    os.remove(file_name)


with DAG(
    dag_id='extract_transform_save_s3',
    default_args={'owner': 'airflow', 'retries': 1},
    description='extract SMILES data, transform and upload to S3',
    schedule_interval='@daily',
    start_date=days_ago(1),
    catchup=False,
) as dag:
    extract_task = PythonOperator(
        task_id='extract_data',
        python_callable=extract_data,
    )

    transform_task = PythonOperator(
        task_id='transform_data',
        python_callable=transform_data,
    )

    save_task = PythonOperator(
        task_id='save_to_s3',
        python_callable=save_to_s3,
    )

    extract_task >> transform_task >> save_task
