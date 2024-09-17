from fastapi import APIRouter, HTTPException, status
from celery.result import AsyncResult
from sqlalchemy import delete
from sqlalchemy.future import select
from src.config import logger
from src.celery_worker import celery
from src.database import async_session_maker
from src.drugs import models, schemas
from src.tasks import add_task_substructure_search


router = APIRouter(prefix="/drugs", tags=["drugs"])


@router.post("", status_code=status.HTTP_201_CREATED)
async def add_drug(drug_data: schemas.DrugCreate) -> schemas.DrugResponse:
    async with async_session_maker() as session:
        async with session.begin():
            new_drug = models.Drug(**drug_data.dict())
            session.add(new_drug)
            await session.flush()
            new_drug_id = new_drug.id
            await session.commit()
            logger.info(f"Added new drug with ID {new_drug_id}")
            return {'id': new_drug_id, **drug_data.dict()}


@router.get("/search")
async def start_substructure_search(mol: str) -> dict:
    task = add_task_substructure_search.delay(mol)
    return {"task_id": task.id, "status": task.status}


@router.get("/search/{task_id}")
async def get_substructure_search_result(task_id: str) -> dict:
    task_result = AsyncResult(task_id, app=celery)
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {"task_id": task_id, "status": "Task completed", "result": task_result.result}
    else:
        return {"task_id": task_id, "status": task_result.state}


@router.get("/{drug_id}")
async def get_drug_full_data(drug_id: int) -> schemas.DrugResponse:
    async with async_session_maker() as session:
        query = select(models.Drug).filter_by(id=drug_id)
        result = await session.execute(query)
        drug_info = result.scalar_one_or_none()

        if not drug_info:
            logger.warning("Drug not found")
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND, detail="Drug not found")
        logger.info(f"Fetched drug with ID {drug_id}")
        return drug_info


@router.get("")
async def get_all_drugs(limit: int = None) -> list[schemas.DrugResponse]:
    async with async_session_maker() as session:
        query = select(models.Drug)
        if limit is not None:
            query = query.limit(limit)
        drugs = await session.execute(query)
        logger.info("Get all drugs")
        return drugs.scalars().all()
    
# def iterator_get_all_drugs(drugs, limit=10):
#     drugs_size = len(drugs)
#     for i in range(0, drugs_size, limit):
#         yield drugs[i: i + limit]


@router.delete("/{drug_id}")
async def delete_drug_by_id(drug_id: int) -> dict:
    async with async_session_maker() as session:
        async with session.begin():
            query = select(models.Drug).filter_by(id=drug_id)
            result = await session.execute(query)
            drug_to_delete = result.scalar_one_or_none()

            if not drug_to_delete:
                logger.warning("Drug not found")
                raise HTTPException(
                    status_code=status.HTTP_404_NOT_FOUND, detail="Drug not found")

            await session.execute(delete(models.Drug).filter_by(id=drug_id))
            await session.commit()
            logger.info(f"Drug with ID {drug_id} deleted")
            return {"message": f"Drug with ID {drug_id} deleted"}


@router.put("/{drug_id}")
async def update_drug(drug_id: int, drug_data: schemas.DrugUpdate) -> schemas.DrugResponse:
    async with async_session_maker() as session:
        async with session.begin():
            query = select(models.Drug).filter_by(id=drug_id)
            result = await session.execute(query)
            drug_to_update = result.scalar_one_or_none()

            if not drug_to_update:
                logger.warning("Drug not found")
                raise HTTPException(status_code=404, detail="Drug not found")

            update_data = drug_data.dict()
            for key, value in update_data.items():
                setattr(drug_to_update, key, value)

            await session.flush()
            await session.commit()
            logger.info(f"Drug with ID {drug_id} updated")
            return drug_to_update
