#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
cmd+=" -r 9:95473178-95475352,9:95476976-95477951,9:95480251-95482499"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
