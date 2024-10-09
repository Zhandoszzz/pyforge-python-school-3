from sqlalchemy import ForeignKey
from sqlalchemy.orm import Mapped, mapped_column
from src.database import Base, int_pk


class Patient(Base):
    id: Mapped[int_pk]
    fname: Mapped[str]
    lname: Mapped[str]
    drug_id: Mapped[int] = mapped_column(ForeignKey("drugs.id"), nullable=False)

    def __str__(self):
        return (
            f"{self.__class__.__name__}(id={self.id}, "
            f"name={self.name!r},"
            f"drug_id={self.drug_id!r})"
        )

    def __repr__(self):
        return str(self)
