#!/bin/bash

N=3         # system size
MAX_OCC=2   # max occupation photon
OMEGA0=1    # photon frequency
H=0.5       # atom splitting
V=0.5       # rydberg-rydberg interactions
# photon matter coupling
G_LIST=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)
G_LIST=(1.5 3.00)
# cavity decay rate
KAPPA=1
# total time
T=20
# time step
DT=0.01
# max bond dimension
MAXDIM=256

MAIN_DIR=./
BIN=$(pwd)/tn_rydberg_in_leaky_cavity  # absolute path: the script cd-s into SAVEDIR

SAVEDIR=${MAIN_DIR}/rydberg_in_leaky_cavity_N${N}_maxocc${MAX_OCC}
if [ ! -d $SAVEDIR ]; then                                                                                  # verifica che la directory da creare non esista già
    if ! mkdir $SAVEDIR ; then                                                                              # verifica che la creazione della directory non sia fallita
    exit -1
    fi
fi 

cd ${SAVEDIR}

for G in ${G_LIST[*]}
do
    INPUT_FILE="input.txt"

    echo "input" > ${INPUT_FILE}
    echo "{" >> ${INPUT_FILE}
    echo "N="${N} >> ${INPUT_FILE}
    echo "max_occ="${MAX_OCC} >> ${INPUT_FILE}
    echo "V="${V} >> ${INPUT_FILE}
    echo "h="${H} >> ${INPUT_FILE}
    echo "g="${G} >> ${INPUT_FILE}
    echo "kappa="${KAPPA} >> ${INPUT_FILE}
    echo "T="${T} >> ${INPUT_FILE}
    echo "dt="${DT} >> ${INPUT_FILE}
    echo "maxDim="${MAXDIM} >> ${INPUT_FILE}
    echo "}" >> ${INPUT_FILE}

    ${BIN} ${INPUT_FILE}
done
