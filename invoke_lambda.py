import boto3
import json
import os


def invoke_lambda():
    session = boto3.Session(
        aws_access_key_id=os.getenv('AWS_ACCESS_KEY_ID'),
        aws_secret_access_key=os.getenv('AWS_SECRET_ACCESS_KEY'),
        region_name='eu-north-1'
    )

    lambda_client = session.client('lambda')

    payload = {
        "names": ["Python", "World", "Anon"]
    }

    response = lambda_client.invoke(
        FunctionName='HelloStudentFunction',
        InvocationType='RequestResponse',
        Payload=json.dumps(payload)
    )

    response_payload = response['Payload'].read()

    print(response_payload)


invoke_lambda()
