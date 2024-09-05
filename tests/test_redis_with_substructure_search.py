import pytest
from src.utils import get_cached_result


@pytest.mark.asyncio
async def test_redis_substructure_search(client):
    await client.get("drugs/search?mol=C")
    cached_result = get_cached_result("search:C")
    assert [] == cached_result