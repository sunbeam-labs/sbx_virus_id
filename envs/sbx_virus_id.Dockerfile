FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_virus_id_env

COPY envs/sbx_virus_id.yml ./

# Install environment
RUN conda env create --file sbx_virus_id.yml --name sbx_virus_id

ENV PATH="/opt/conda/envs/sbx_virus_id/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_virus_id", "/bin/bash", "-c"]

# Run
CMD "bash"