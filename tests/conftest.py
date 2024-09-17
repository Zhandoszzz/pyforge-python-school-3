from httpx import AsyncClient, ASGITransport
from src.main import app
import pytest_asyncio


@pytest_asyncio.fixture
async def client():
    async with AsyncClient(transport=ASGITransport(app=app), base_url="http://test") as client:
        yield client
