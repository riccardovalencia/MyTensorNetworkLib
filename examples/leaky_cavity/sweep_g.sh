#!/bin/bash
# Run leaky_cavity for several couplings g (units of the critical coupling): each run writes into its
# own folder data/<run>/ (the run name contains g), plot it with plot_leaky_cavity.py <run>.
# Usage: ./sweep_g.sh   (from examples/leaky_cavity, after building the examples)

N=3          # number of sites (cavity + N-1 spins)
MAX_OCC=2    # maximal photon occupation
H=0.5        # atomic splitting
V=0.5        # Rydberg interaction
KAPPA=1      # cavity decay rate
T=20         # total time
DT=0.01      # time step
MAX_DIM=256  # maximal bond dimension
G_LIST=(1.5 3.0)

INPUT_FILE=$(mktemp)
trap 'rm -f "$INPUT_FILE"' EXIT

for G in "${G_LIST[@]}"
do
    cat > "$INPUT_FILE" <<EOF
input
{
N = $N
max_occ = $MAX_OCC
h = $H
V = $V
g = $G
kappa = $KAPPA
T = $T
dt = $DT
max_dim = $MAX_DIM
}
EOF
    ./leaky_cavity "$INPUT_FILE" || exit 1
done
