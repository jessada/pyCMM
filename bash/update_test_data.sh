#!/bin/bash
source $PYCMM/bash/cmm_functions.sh

test_data_root="$PYCMM/pycmm"

cmd="$test_data_root/proc/db/test/data/connector/SQLiteDB/test_anno_col_to_db_col_2/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/db/test/data/connector/SQLiteDB/test_sample_id_to_db_col_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/db/test/data/connector/SQLiteDB/test_sample_id_to_tbl_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/SQLiteDBReader/test_filter_non_intergenic_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/QryCall/test_shared_mutation_2/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/QryCall/test_shared_mutation_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/QryRecord/test_pathogenic_count_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/QryRecord/test_parse_intervar_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/QryRecord/test_parse_exac_constraint_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/QryRecord/test_cal_est_kvot_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader/QryRecord/test_max_ref_maf_2/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_coloring_zygosity_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_color_genes_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_select_dataset_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_color_genes_2/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_filter_genes_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_show_shared_variants_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_coloring_zygosity_2/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_header_corrections_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_filter_non_downstream_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_filter_non_intergenic_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_filter_non_intronic_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_filter_non_synonymous_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_filter_non_upstream_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_filter_non_utr_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_max_ref_maf_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_parse_exac_constraint_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_parse_intervar_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/dbreader_xls/QryRecordXls/test_pathogenic_count_xls_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_unicode_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_multiple_regions_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_datasets_1/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_datasets_2/gen_data.sh"
eval_cmd "$cmd"

cmd="$test_data_root/proc/mutrep/test/data/mutrep/MutRepController/test_filter_non_recessive_gene_1/gen_data.sh"
eval_cmd "$cmd"
