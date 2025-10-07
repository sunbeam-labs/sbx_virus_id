FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_cenote_taker_env

COPY envs/sbx_cenote_taker.yml ./

# Install environment
RUN conda env create --file sbx_cenote_taker.yml --name sbx_cenote_taker

ENV PATH="/opt/conda/envs/sbx_cenote_taker/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_cenote_taker", "/bin/bash", "-c"]

# Run
CMD "bash"
