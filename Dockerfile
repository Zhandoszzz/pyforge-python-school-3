FROM continuumio/miniconda3

WORKDIR /src

RUN conda install -c conda-forge rdkit \
    && pip install fastapi uvicorn

COPY ./src /src

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]