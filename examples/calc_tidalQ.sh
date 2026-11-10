#!/bin/sh
set -e

script_dir=$(CDPATH= cd -P -- "$(dirname -- "$0")" && pwd)
solver=${TIDALQ:-"$script_dir/../calc_tidalQ"}
model_dir=${1:-"$script_dir"}
output_dir=${2:-sample-results}

mkdir -p -- "$output_dir"

for m2_model in "$model_dir"/model_tidal_Q_model*_1.dat; do
    id=${m2_model##*/model_tidal_Q_model}
    id=${id%_1.dat}
    nodal_model=${m2_model%_1.dat}_2.dat

    m2=$("$solver" -W1.4056e-4 -M2 -S"$m2_model")
    nodal=$("$solver" -W1.065e-8 -M0 -S"$nodal_model")
    printf '%s\n%s\n' "$m2" "$nodal" > "$output_dir/tidalQ.$id.out"
done
