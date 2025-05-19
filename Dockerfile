# first stage
FROM python:3.13 AS builder

COPY requirements.txt .
COPY requirements/* ./requirements/

# install dependencies to the local user directory (eg. /root/.local)
RUN pip install --user -r requirements.txt

# second unnamed stage
FROM python:3.13-slim
WORKDIR /code

# copy only the dependencies installation from the 1st stage image
COPY --from=builder /root/.local/bin /root/.local
COPY . .

# update PATH environment variable
ENV PATH=/root/.local:$PATH

# install the code
RUN pip install -e . -r requirements.txt
RUN chmod -R +rx /code

ENTRYPOINT [ "fastq_demux" ]
CMD [ "--help" ]
