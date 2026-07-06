#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$SCRIPT_DIR/build/duneanaobj/lib:$SCRIPT_DIR/build/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$SCRIPT_DIR/duneanaobj:$ROOT_INCLUDE_PATH

if [ ! "$#" -ge 1 ]; then
    echo "Usage: ./make_caf_simulation.sh [input_files]"
fi

offset=0
for i in "$@"
do
    output=`sed 's/.\{3\}$//' <<< "${i}"`.caf.root
    "$(dirname "$0")/build/make_standalone_simulation" $output $offset $i
    ((offset+=1000))
done