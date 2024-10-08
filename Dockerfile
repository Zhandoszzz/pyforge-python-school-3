FROM continuumio/miniconda3

WORKDIR /src

COPY ./requirements.txt .

RUN pip install --no-cache-dir -r /src/requirements.txt

COPY ./src /src

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]