#!/bin/bash
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-r {name}          query region(s) (required)
-o {file}          output direcotry (required)
EOF
)

while getopts ":r:o:" OPTION; do
  case "$OPTION" in
    r)
      qry_regions_txt="$OPTARG"
      ;;
    o)
      out_dir="$OPTARG"
      ;;
    *)
      die "unrecognized option from executing: $0 $@"
      ;;
  esac
done

[ ! -z $qry_regions_txt ] || die "Please query regions (-r)"
[ ! -z $out_dir ] || die "Please specify output folder (-o)"
[ -d $out_dir ] || die "$out_dir is not a valid folder name"

tmp_dir="$HOME/tmp_db"
mkdir -p "$tmp_dir"
## ****************************************  display configuration  ****************************************
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "parameters"
info_msg "  $params"
info_msg
info_msg "description"
info_msg "  A script to generate sqliteDB test data"
info_msg
info_msg
## display required configuration
info_msg "overall configuration"
display_param "region(s) (-r)"
IFS=$',' read -ra qry_regions_list <<< "$qry_regions_txt"
if [ $((${#qry_regions_list[@]})) -gt 1 ]; then
for (( i=0; i<$((${#qry_regions_list[@]})); i++ ))
do
    display_param "  region $(( i+1 ))" "${qry_regions_list[$i]}"
done
fi
display_param "output folder (-o)" "$out_dir"
display_param "tmp folder" "$tmp_dir"

# ****************************************  executing  ****************************************

new_section_txt "E X E C U T I N G"

source_db="$SQLITE_DB_GRCH37"

schema_sql="$tmp_dir/schema.sql"
dump_sql="$tmp_dir/dump.sql"
insert_all_gtz_annos_sql="$tmp_dir/insert_all_gtz_annos.sql"
insert_cmm_tables_sql="$tmp_dir/insert_cmm_tables.sql"

new_db="$out_dir/input.db"

# dump schema 
cmd="sqlite3 $source_db .schema > $schema_sql"
eval_cmd "$cmd"

> $dump_sql
for (( n=0; n<$((${#qry_regions_list[@]})); n++ ))
do
    IFS=$':' read -ra region_split <<< "${qry_regions_list[$n]}"
    chrom="${region_split[0]}"
    IFS=$'-' read -ra pos_split <<< "${region_split[1]}"
    start_pos="${pos_split[0]}"
    end_pos="${pos_split[1]}"

    # dump all_gtz_annos table by region
    cmd="echo .mode insert all_gtz_annos >> $dump_sql"
    eval_cmd "$cmd"
    cmd="echo \".output $insert_all_gtz_annos_sql$n\" >> $dump_sql"
    eval_cmd "$cmd"
    cmd="echo \"select * from all_gtz_annos where CHROM = '$chrom' and POS >= $start_pos and POS <= $end_pos;\" >> $dump_sql"
    eval_cmd "$cmd"
done

# dump cmm_tables table
cmd="echo .mode insert cmm_tables >> $dump_sql"
eval_cmd "$cmd"
cmd="echo \".output $insert_cmm_tables_sql\" >> $dump_sql"
eval_cmd "$cmd"
cmd="echo \"select * from cmm_tables;\" >> $dump_sql"
eval_cmd "$cmd"
cmd="sqlite3 $source_db < $dump_sql"
eval_cmd "$cmd"

#create new database
> $new_db
cmd="sqlite3 $new_db < $schema_sql"
eval_cmd "$cmd"
cmd="sqlite3 $new_db < $insert_cmm_tables_sql"
eval_cmd "$cmd"
for (( n=0; n<$((${#qry_regions_list[@]})); n++ ))
do
    cmd="sqlite3 $new_db < $insert_all_gtz_annos_sql$n"
    eval_cmd "$cmd"
done

#rm -r $tmp_dir

new_section_txt "F I N I S H <$script_name>"
