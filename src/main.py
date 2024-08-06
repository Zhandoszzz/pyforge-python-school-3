from os import getenv
from pydantic import BaseModel
from rdkit import Chem
from fastapi import FastAPI, HTTPException, status, File, UploadFile
import json


app = FastAPI()
molecules_db = {}


class Molecule(BaseModel):
    id: int
    smile: str


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID")}


@app.post("/molecules", status_code=status.HTTP_201_CREATED)
def add_molecule(mol: Molecule):
    if mol.id in molecules_db:
        raise HTTPException(status_code=status.HTTP_409_CONFLICT,
                            detail='Molecule already exists')
    molecules_db[mol.id] = mol.smile
    return mol


@app.get("/molecules/{mol_id}")
def get_molecule_by_id(mol_id: int):
    if mol_id in molecules_db:
        return molecules_db[mol_id]
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                        detail='Molecule Not Found')


@app.put("/molecules/{mol_id}")
def update_molecule_by_id(mol_id: int, updated_mol: Molecule):
    if mol_id in molecules_db:
        del molecules_db[mol_id]
        molecules_db[updated_mol.id] = updated_mol.smile
        return updated_mol
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                        detail='Molecule Not Found')


@app.delete("/molecules/{mol_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_molecule_by_id(mol_id: int):
    if mol_id in molecules_db:
        del molecules_db[mol_id]
        return f"Molecule with id {mol_id} is deleted"
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                        detail='Molecule Not Found')


@app.get("/molecules")
def list_all_molecules():
    return molecules_db


@app.get("/search")
def substructure_search(mol: str):
    mol = Chem.MolFromSmiles(mol)
    result = {id: smile for id, smile in molecules_db.items() if Chem.MolFromSmiles(
        smile).HasSubstructMatch(mol)}
    return result


@app.post("/upload", status_code=status.HTTP_201_CREATED)
async def upload_file_with_molecules(file: UploadFile = File(...)):
    if file.content_type != 'text/plain':
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST,
                            detail='Expected txt file')

    content = await file.read()
    content_str = content.decode('utf-8')
    try:
        molecules_list = json.loads(content_str)
        for molecule in molecules_list:
            molecule_id = molecule['id']
            molecule_smile = molecule['smile']
            if molecule_id in molecules_db:
                raise HTTPException(status_code=status.HTTP_409_CONFLICT,
                                    detail='Molecule already exists')
            molecules_db[molecule_id] = molecule_smile
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST,
                            detail='Error processing file')

    return molecules_db
