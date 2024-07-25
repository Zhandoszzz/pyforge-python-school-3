from pydantic import BaseModel
from rdkit import Chem
from fastapi import FastAPI, HTTPException, status, File, UploadFile
import json


app = FastAPI()
molecules_db = []


class Molecule(BaseModel):
    id: int
    smile: str


@app.post("/molecules", status_code=status.HTTP_201_CREATED)
def add_molecule(mol: Molecule):
    molecules_db.append(mol.dict())
    return mol


@app.get("/molecules/{mol_id}")
def get_molecule_by_id(mol_id: int):
    for mol in molecules_db:
        if mol['id'] == mol_id:
            return mol
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                        detail='Molecule Not Found')


@app.put("/molecules/{mol_id}")
def update_molecule_by_id(mol_id: int, updated_mol: Molecule):
    for index, mol in enumerate(molecules_db):
        if mol['id'] == mol_id:
            molecules_db[index] = updated_mol.dict()
            return updated_mol
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                        detail='Molecule Not Found')


@app.delete("/molecules/{mol_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_molecule_by_id(mol_id: int):
    for index, mol in enumerate(molecules_db):
        if mol['id'] == mol_id:
            molecules_db.pop(index)
            return f"Molecule with id {mol_id} is deleted"
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                        detail='Molecule Not Found')


@app.get("/molecules")
def list_all_molecules():
    return molecules_db


@app.get("/search")
def substructure_search(mol: str):
    result = []
    mol = Chem.MolFromSmiles(mol)
    result = [m['smile'] for m in molecules_db if Chem.MolFromSmiles(
        m['smile']).HasSubstructMatch(mol)]
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
        molecules_db.extend(molecules_list)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST,
                            detail='Error processing file')

    return molecules_db
