#!/bin/bash

# ==============================================================
# === Declarations
# Make sure script exits if error encountered
set -o pipefail
set -o errexit

# Set defaults for options
RMPNG=0
RMDAT=0
FRAMERATE=3


print_help() {
    echo "Cleans, makes, and runs simulation. Then, it creates animations."
    echo
    echo "-h     : Prints this message."
    echo "-r     : (default=false) to clean directory of png files"
    echo "-f OPT : (default=3) sets the frames per second of output video" 
    echo "-c     : (default=false) to clean data directory of *.csv files"
}

while getopts o:f:rh flag
do
    case "${flag}" in
        c) RMDAT=1;;
        f) FRAMERATE="${OPTARG}";;
        r) RMPNG=1;;
        h) print_help
            exit 1;;
    esac
done

# ==============================================================
# === Main script
# Make and run simulation.
make clean && make && echo $'\nStarting program..\n\n' && ./fdtd

# Plot outputs and create mp4 files.
echo $'\n\nStarting plotting..'
(python ./plots.py FRAMERATE RMPNG)


if [[ $RMDAT = 'true' ]]
then
    echo $'Cleaning data..'
    rm data/*.csv
fi

echo $'\nExiting\n'
