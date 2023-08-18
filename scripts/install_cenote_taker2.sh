#!/bin/bash

cd {params.ext_fp}
if [ ! -d "Cenote-Taker2" ]; then
    git clone https://github.com/mtisza1/Cenote-Taker2.git
else
    echo "Cenote-Taker2 directory already exists"
fi

#if {{ conda env list | grep 'cenote-taker2_env'; }} >/dev/null 2>&1; then
#    echo "Cenote-Taker2 env already exists"
#else
#    conda create --name cenote-taker2_env --file cenote_taker2_env_explicit.txt
#fi

#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh
#conda activate cenote-taker2_env

pip install phanotate

# check if all dbs are installed already


# with all the options (75GB). The PDB database (--hhPDB) takes about 2 hours to download.
if [ -z {params.db_fp} ]; then
    cd Cenote-Taker2/
    if [ -d "pdb70/" ]; then
        echo "CenoteTaker2 DBs already installed"
    else
        python update_ct2_databases.py --hmm True --protein True --rps True --taxdump True --hhCDD True --hhPFAM True --hhPDB True >> {log}
    fi
else
    mkdir -p {params.db_fp}
    cd {params.db_fp}
    python {params.ext_fp}/Cenote-Taker2/update_ct2_databases.py --hmm True --protein True --rps True --taxdump True --hhCDD True --hhPFAM True --hhPDB True >> {log}
fi

touch {output}