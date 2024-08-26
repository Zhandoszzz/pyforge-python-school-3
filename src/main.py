from os import getenv
from fastapi import FastAPI
from src.drugs.router import router as drug_router
from src.patients.router import router as patient_router


app = FastAPI()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID")}


app.include_router(drug_router)
app.include_router(patient_router)
