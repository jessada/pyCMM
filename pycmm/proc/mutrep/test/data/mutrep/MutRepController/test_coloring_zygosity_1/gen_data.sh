#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
regions="9:95473178-95477951"
regions+=",9:95480251-95482499"
cmd+=" -r $regions"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
