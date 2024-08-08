import pytest
from fastapi.testclient import TestClient
from src.main import app


def add_molecules(client):
    molecules = [
        {"id": 1, "smile": "CCO"},
        {"id": 2, "smile": "c1ccccc1"},
        {"id": 3, "smile": "CC(=O)O"},
        {"id": 4, "smile": "CC(=O)Oc1ccccc1C(=O)O"}
    ]
    for mol in molecules:
        client.post("/molecules", json=mol)


@pytest.fixture(scope="module")
def client():
    client = TestClient(app)
    add_molecules(client)
    return client


@pytest.mark.parametrize(
    "query_mol, expected_result",
    [
        ("c1ccccc1", {"2": "c1ccccc1", "4": "CC(=O)Oc1ccccc1C(=O)O"}),
        ("CC", {"1": "CCO", "3": "CC(=O)O", "4": "CC(=O)Oc1ccccc1C(=O)O"}),
        ("CC=O", {"3": "CC(=O)O", "4": "CC(=O)Oc1ccccc1C(=O)O"})
    ]
)
def test_substructure_search_valid(client, query_mol, expected_result):
    response = client.get(f'/search?mol={query_mol}')
    assert response.status_code == 200
    assert response.json() == expected_result


def test_substructure_search_no_match(client):
    response = client.get('/search?mol=NCCN')
    assert response.status_code == 200
    assert response.json() == {}


def test_substructure_search_invalid_mol(client):
    response = client.get('/search?mol=hello')
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid input"}
