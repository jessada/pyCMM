#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
regions="3:133370982-133371807"
regions+=",3:133373894-133375034"
cmd+=" -r $regions"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
