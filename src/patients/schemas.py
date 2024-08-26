from pydantic import BaseModel


class PatientBase(BaseModel):
    fname: str
    lname: str
    drug_id: int


class PatientCreate(PatientBase):
    pass


class Patient(PatientBase):
    id: int

    class Config:
        orm_mode = True