from pydantic import BaseModel


class DrugBase(BaseModel):
    name: str
    smiles: str


class DrugCreate(DrugBase):
    pass


class DrugUpdate(DrugBase):
    pass


class Drug(DrugBase):
    id: int

    class Config:
        orm_mode = True
