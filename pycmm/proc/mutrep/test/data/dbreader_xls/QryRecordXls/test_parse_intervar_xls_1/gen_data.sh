#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
cmd+=" -r 1:19018405-19018405,1:36939403-36939403,1:43812237-43812237,1:43812255-43812255,2:38302361-38302361,3:133371424-133371424"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
