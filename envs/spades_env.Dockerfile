FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_virus_id_env

COPY envs/spades_env.yml ./

# Install environment
RUN conda env create --file spades_env.yml --name spades

ENV PATH="/opt/conda/envs/spades/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "spades", "/bin/bash", "-c"]

# Run
CMD "bash"