FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_virus_id_env

COPY envs/cenote_taker_env.yml ./

# Install environment
RUN conda env create --file cenote_taker_env.yml --name cenote_taker

ENV PATH="/opt/conda/envs/cenote_taker/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "cenote_taker", "/bin/bash", "-c"]

# Run
CMD "bash"