# Ensure ups is available if not already sourced
if ! command -v setup &> /dev/null; then
    echo "Warning: 'setup' command not found. You may need to source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh first."
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export HDF5_HOME=$SCRIPT_DIR/hdf5_lib/

setup cmake v3_21_4
setup duneanaobj v03_11_00 -q e26:prof
