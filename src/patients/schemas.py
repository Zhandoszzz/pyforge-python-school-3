from pydantic import BaseModel


class PatientBase(BaseModel):
    fname: str
    lname: str
    drug_id: int


class PatientCreate(PatientBase):
    pass


class PatientResponse(PatientBase):
    id: int

    class Config:
        from_attributes = True