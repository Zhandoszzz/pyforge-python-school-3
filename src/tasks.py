from src.celery_worker import celery
from src.database import async_session_maker
from src.utils import set_cache, get_cached_result
from src.config import logger
from fastapi import HTTPException, status
from sqlalchemy.future import select
from src.drugs import models
from rdkit import Chem


@celery.task
async def add_task_substructure_search(mol):
    cache_key = f"search:{mol}"
    cached_result = get_cached_result(cache_key)

    if cached_result is not None:
        logger.info("Return cached result")
        return cached_result

    async with async_session_maker() as session:
        mol = Chem.MolFromSmiles(mol)
        query = select(models.Drug)
        query_result = await session.execute(query)
        drugs = query_result.scalars().all()
        if not mol:
            logger.warning("Invalid SMILE")
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST,
                                detail='Invalid input')
        result = [drug for drug in drugs if Chem.MolFromSmiles(
            drug.smiles).HasSubstructMatch(mol)]
        set_cache(cache_key, result)
        logger.info(f"Set cache - {cache_key}")
        logger.info("Performed substructure search")
        return result