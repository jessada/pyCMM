#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
regions="2:179408713-179411011"
regions+=",6:31324206-31324601"
regions+=",6:31740760-31740763"
regions+=",7:100550486-100550486"
regions+=",19:17392630-17392630"
regions+=",19:23159486-23159486"
regions+=",19:32083223-32083250"
cmd+=" -r $regions"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
