FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_virus_id_env

COPY envs/virsorter_env.yml ./

# Install environment
RUN conda env create --file virsorter_env.yml --name virsorter

ENV PATH="/opt/conda/envs/virsorter/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "virsorter", "/bin/bash", "-c"]

# Run
CMD "bash"