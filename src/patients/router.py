from fastapi import APIRouter, Depends, HTTPException, status
from src.database import async_session_maker
from sqlalchemy.future import select
from sqlalchemy import delete
from . import models, schemas

router = APIRouter(prefix="/patients", tags=["patients"])


@router.post("", status_code=status.HTTP_201_CREATED)
async def add_patient(patient_data: schemas.PatientCreate):
    async with async_session_maker() as session:
        async with session.begin():
            new_patient = models.Patient(**patient_data.dict())
            session.add(new_patient)
            await session.flush()
            new_patient_id = new_patient.id
            await session.commit()
            return {'id': new_patient_id, **patient_data.dict()}


@router.delete("/{patient_id}")
async def delete_patient_by_id(patient_id: int):
    async with async_session_maker() as session:
        async with session.begin():
            query = select(models.Patient).filter_by(id=patient_id)
            result = await session.execute(query)
            patient_to_delete = result.scalar_one_or_none()

            if not patient_to_delete:
                raise HTTPException(
                    status_code=status.HTTP_404_NOT_FOUND, detail="Patient not found")

            await session.execute(delete(models.Patient).filter_by(id=patient_id))
            await session.commit()
            return {"message": f"Patient with ID {patient_id} deleted"}


@router.get("")
async def get_all_patients():
    async with async_session_maker() as session:
        query = select(models.Patient)
        drugs = await session.execute(query)
        return drugs.scalars().all()
