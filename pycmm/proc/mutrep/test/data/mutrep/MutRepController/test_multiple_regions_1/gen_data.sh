#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
regions="10:89601588-89601588"
regions+=",10:89617802-89618430"
regions+=",10:89682209-89682833"
regions+=",11:57427955-57427955"
cmd+=" -r $regions"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
