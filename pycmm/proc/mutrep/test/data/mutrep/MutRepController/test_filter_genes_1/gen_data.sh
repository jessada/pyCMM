#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
regions="9:95468981-95474259"
regions+=",9:95477131-95477873"
regions+=",9:95526955-95529019"
regions+=",9:95570290-95571939"
cmd+=" -r $regions"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
