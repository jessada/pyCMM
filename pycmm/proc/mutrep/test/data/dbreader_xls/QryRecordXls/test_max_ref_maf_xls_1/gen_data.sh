#!/bin/bash

gen_test_db_script="$PYCMM/bash/gen_test_db.sh"

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cmd="$gen_test_db_script"
regions="1:1268000-1268000"
regions+=",1:1334052-1334052"
regions+=",1:1458900-1458900"
regions+=",1:1857306-1857306"
regions+=",1:3669356-3669356"
regions+=",1:3683159-3683159"
regions+=",1:47904668-47904668"
regions+=",12:21176135-21176135"
regions+=",17:16635980-16635980"
regions+=",17:16675944-16675944"
regions+=",17:17039561-17039561"
cmd+=" -r $regions"
cmd+=" -o $script_dir"
eval_cmd "$cmd"
