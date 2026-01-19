FROM python:3.11-slim

WORKDIR /app

RUN pip install --no-cache-dir numpy matplotlib

COPY MoC.py .

ENV MPLBACKEND=Agg

CMD ["python", "MoC.py"]