#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
cmd+=" -r 6:78171941-78172992"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
