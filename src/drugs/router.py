from fastapi import APIRouter, HTTPException, status
from src.database import async_session_maker
from sqlalchemy.future import select
from sqlalchemy import delete
from rdkit import Chem
from . import models, schemas


router = APIRouter(prefix="/drugs", tags=["drugs"])


@router.post("", status_code=status.HTTP_201_CREATED)
async def add_drug(drug_data: schemas.DrugCreate):
    async with async_session_maker() as session:
        async with session.begin():
            new_drug = models.Drug(**drug_data.dict())
            session.add(new_drug)
            await session.flush()
            new_drug_id = new_drug.id
            await session.commit()
            return {'id': new_drug_id, **drug_data.dict()}


@router.get("/search")
async def substructure_search(mol: str):
    async with async_session_maker() as session:
        mol = Chem.MolFromSmiles(mol)
        query = select(models.Drug)
        query_result = await session.execute(query)
        drugs = query_result.scalars().all()
        if not mol:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST,
                                detail='Invalid input')
        result = [drug for drug in drugs if Chem.MolFromSmiles(
            drug.smiles).HasSubstructMatch(mol)]
        return result


@router.get("/{drug_id}")
async def get_drug_full_data(drug_id: int):
    async with async_session_maker() as session:
        query = select(models.Drug).filter_by(id=drug_id)
        result = await session.execute(query)
        drug_info = result.scalar_one_or_none()

        if not drug_info:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND, detail="Drug not found")

        return drug_info


@router.get("")
async def get_all_drugs():
    async with async_session_maker() as session:
        query = select(models.Drug)
        drugs = await session.execute(query)
        return drugs.scalars().all()


@router.delete("/{drug_id}")
async def delete_drug_by_id(drug_id: int):
    async with async_session_maker() as session:
        async with session.begin():
            query = select(models.Drug).filter_by(id=drug_id)
            result = await session.execute(query)
            drug_to_delete = result.scalar_one_or_none()

            if not drug_to_delete:
                raise HTTPException(
                    status_code=status.HTTP_404_NOT_FOUND, detail="Drug not found")

            await session.execute(delete(models.Drug).filter_by(id=drug_id))
            await session.commit()
            return {"message": f"Drug with ID {drug_id} deleted"}


@router.put("/{drug_id}")
async def update_drug(drug_id: int, drug_data: schemas.DrugUpdate):
    async with async_session_maker() as session:
        async with session.begin():
            query = select(models.Drug).filter_by(id=drug_id)
            result = await session.execute(query)
            drug_to_update = result.scalar_one_or_none()

            if not drug_to_update:
                raise HTTPException(status_code=404, detail="Drug not found")

            update_data = drug_data.dict()
            for key, value in update_data.items():
                setattr(drug_to_update, key, value)

            await session.flush()
            await session.commit()

            return drug_to_update
