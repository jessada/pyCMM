from pycmm.utils import mylogger
import pkg_resources
import sys

param_display_fmt = "  {name:<45}{value}"

def new_section_txt(txt):
    adj_txt = " " + txt + " "
    mylogger.info("")
    mylogger.info(adj_txt.center(140,"*"))

def disp_header(header_txt):
    mylogger.info("")
    mylogger.info(header_txt)

def disp_subheader(subheader_txt):
    mylogger.info("  " + subheader_txt)

def disp_param(param_name, param_value):
    mylogger.info(param_display_fmt.format(name=param_name+":",
                                           value=param_value))

def debug_param(param_name, param_value):
    mylogger.debug(param_display_fmt.format(name=param_name+":",
                                            value=param_value))

def disp_subparam(subparam_name, subparam_value):
    disp_param("  "+subparam_name, subparam_value)

def show_config(app_description,
                time_stamp,
                required_params,
                optional_params,
                ):
    disp_header("Description")
    mylogger.info("  " + app_description)
    disp_header("Version and environment configuration")
    disp_param("pyCMM version", pkg_resources.get_distribution("pycmm").version)
    disp_param("parameters", " ".join(sys.argv[1:]))
    debug_param("debug mode", "ON")
    disp_param("time stamp", time_stamp)
    disp_header("Required parameters")
    for key in required_params:
        disp_param(key, required_params[key])
    disp_header("Optional parameters")
    for key in optional_params:
        disp_param(key, optional_params[key])
### ****************************************  display configuration  ****************************************
#new_section_txt(" S T A R T <" + script_name + "> ")
#mylogger.info("")
#disp_header("script configuration")
#disp_param("parameters", " ".join(sys.argv[1:]))
#disp_param("timestamp", datetime.datetime.now())
#mylogger.info("")
#
### display required configuration
#disp_header("required configuration")
#disp_param("xls output file (-o)", out_file)
#disp_param("master columns count (-o)", n_master_cols)
#mylogger.info("")
#
### display csvs configuration
#disp_header("csvs configuration (-s)(" + str(len(csvs_list)) + " sheet(s))")
#for i in xrange(len(csvs_list)):
#    (sheet_name, sheet_csv) = csvs_list[i].split(',')
#    disp_param("sheet name #"+str(i+1), sheet_name)
#    disp_param("sheet csv  #"+str(i+1), sheet_csv)
#mylogger.info("")
#
### display optional configuration
#disp_header("optional configuration")
#if len(addn_csvs_list) > 0:
#    n_addn_sheet = len(addn_csvs_list)
#    header_txt = "additional csv sheets configuration (-A)"
#    header_txt += "(" + str(n_addn_sheet) + " sheet(s))"
#    disp_subheader(header_txt)
#    for i in xrange(len(addn_csvs_list)):
#        (sheet_name, sheet_csv) = addn_csvs_list[i].split(',')
#        disp_subparam("sheet name #"+str(i+1), sheet_name)
#        disp_subparam("sheet csv  #"+str(i+1), sheet_csv)
#if marked_key_range is not None :
#    disp_subheader("marked key range")
#    disp_subparam("start key", marked_start_key)
#    disp_subparam("end key", marked_end_key)
#if len(frequency_ratios) > 0:
#    disp_subheader("frequency_ratios (-F)")
#    for i in xrange(len(frequency_ratios)):
#	(col_name, freq) = frequency_ratios[i].split(':')
#        disp_subparam(col_name, freq)
#if len(inc_criteria) > 0:
#    disp_param("inclusion criteria", ",".join(inc_criteria))
#if len(xtra_attribs) > 0:
#    disp_subheader("extra attributes (-E)")
#    for i in xrange(len(xtra_attribs)):
#        disp_subparam("extra attributes #"+str(i+1), xtra_attribs[i])
#disp_subheader("zygosity codes (-Z)")
#for zygo_key in ZYGO_CODES:
#    disp_subparam(zygo_key, ZYGO_CODES[zygo_key])
#if args.cell_colors is not None:
#    disp_param("cell colors (-K)", args.cell_colors)
#if len(color_region_infos) > 0:
#    disp_subheader("color regions information (-C)")
#    for i in xrange(len(color_region_infos)):
#        color_region_info = color_region_infos[i]
#        disp_subparam("color info #"+str(i+1), color_region_info.raw_info)
#if dev_mode:
#    disp_param("developer mode (-D)", "ON")
#
#
### ****************************************  executing  ****************************************
#def set_layout(ws, record_size, col_idx_mg):
#    # hide key, end postion and effect predictors columns
#    ws.set_column(col_idx_mg.IDX_KEY, col_idx_mg.IDX_KEY, None, None, {'hidden': True})
#    ws.set_column(col_idx_mg.IDX_END, col_idx_mg.IDX_END, None, None, {'hidden': True})
##    ws.set_column(col_idx_mg.IDX_PL, col_idx_mg.IDX_MTPRED, None, None, {'hidden': True})
#    # set column width
#    ws.set_column(col_idx_mg.IDX_FUNC, col_idx_mg.IDX_FUNC, 12)
#    ws.set_column(col_idx_mg.IDX_GENE, col_idx_mg.IDX_GENE, 10)
#    ws.set_column(col_idx_mg.IDX_OAF, col_idx_mg.IDX_1000G, 7)
#    ws.set_column(col_idx_mg.IDX_DBSNP, col_idx_mg.IDX_DBSNP, 10)
#    ws.set_column(col_idx_mg.IDX_CHR, col_idx_mg.IDX_CHR, 2)
#    ws.set_column(col_idx_mg.IDX_START, col_idx_mg.IDX_START, 10)
#    ws.set_column(col_idx_mg.IDX_REF, col_idx_mg.IDX_OBS, 6)
#    # freeze panes
#    ws.freeze_panes(HORIZONTAL_SPLIT_IDX, col_idx_mg.IDX_PL)
#    # set auto filter
#    ws.autofilter(0, 0, 0, record_size-1)
#
#def write_header(ws,
#                 cell_fmt_mg,
#                 header_rec,
#                 rec_size,
#                 col_idx_mg,
#                 xtra_attribs,
#                 ):
#    cell_fmt = cell_fmt_mg.cell_fmts[DFLT_FMT]
#    for col_idx in xrange(rec_size):
#        ws.write(0, col_idx, header_rec[col_idx], cell_fmt)
#    ws.write(0, col_idx_mg.IDX_1000G, '1000G', cell_fmt)
#    ws.write(0, col_idx_mg.IDX_ESP6500, 'ESP6500', cell_fmt)
#    ws.write(0, col_idx_mg.IDX_DBSNP, 'dbSNP', cell_fmt)
#    ws.write(0, col_idx_mg.IDX_START, 'start position', cell_fmt)
#    ws.write(0, col_idx_mg.IDX_END, 'end position', cell_fmt)
#    for attrib_idx in xrange(len(xtra_attribs)):
#        ws.write(0, rec_size + attrib_idx, xtra_attribs[attrib_idx], cell_fmt)
#
#def write_content(ws,
#                  cell_fmt_mg,
#                  row,
#                  content_rec,
#                  rec_size,
#                  col_idx_mg,
#                  xtra_attribs,
#                  ):
#    dflt_cell_fmt = cell_fmt_mg.cell_fmts[DFLT_FMT]
#    rare = content_rec.is_rare
#    if rare:
#        cell_fmt = cell_fmt_mg.cell_fmts['YELLOW']
#    else:
#        cell_fmt = dflt_cell_fmt
#    marked_color = content_rec.marked_color
#    if marked_color is not None:
#        marked_fmt = cell_fmt_mg.cell_fmts[marked_color]
#    else:
#        marked_fmt = cell_fmt
#    ws.write(row, col_idx_mg.IDX_KEY, content_rec.key, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_FUNC, content_rec.func, marked_fmt)
#    ws.write(row, col_idx_mg.IDX_GENE, content_rec.gene, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_EXFUNC, content_rec.ex_func, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_AACHANGE, content_rec.aa_change, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_OAF, str(content_rec.oaf), cell_fmt)
#    ws.write(row, col_idx_mg.IDX_1000G, str(content_rec.maf), cell_fmt)
#    ws.write(row, col_idx_mg.IDX_ESP6500, str(content_rec.esp6500), cell_fmt)
#    ws.write(row, col_idx_mg.IDX_DBSNP, content_rec.dbsnp, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_CHR, content_rec.chrom, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_START, content_rec.start, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_END, content_rec.end, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_REF, content_rec.ref, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_OBS, content_rec.obs, cell_fmt)
#    ws.write(row, col_idx_mg.IDX_PL, content_rec.pl, dflt_cell_fmt)
#    if content_rec.pl_harmful:
#        harmful_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_HARMFUL]]
#    else:
#        harmful_fmt = dflt_cell_fmt
#    ws.write(row, col_idx_mg.IDX_PLPRED, content_rec.pl_pred, harmful_fmt)
#    ws.write(row, col_idx_mg.IDX_SIFT, content_rec.sift, dflt_cell_fmt)
#    if content_rec.sift_harmful:
#        harmful_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_HARMFUL]]
#    else:
#        harmful_fmt = dflt_cell_fmt
#    ws.write(row, col_idx_mg.IDX_SIFTPRED, content_rec.sift_pred, harmful_fmt)
#    ws.write(row, col_idx_mg.IDX_PP, content_rec.pp, dflt_cell_fmt)
#    if content_rec.pp_harmful:
#        harmful_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_HARMFUL]]
#    else:
#        harmful_fmt = dflt_cell_fmt
#    ws.write(row, col_idx_mg.IDX_PPPRED, content_rec.pp_pred, harmful_fmt)
#    ws.write(row, col_idx_mg.IDX_LRT, content_rec.lrt, dflt_cell_fmt)
#    if content_rec.lrt_harmful:
#        harmful_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_HARMFUL]]
#    else:
#        harmful_fmt = dflt_cell_fmt
#    ws.write(row, col_idx_mg.IDX_LRTPRED, content_rec.lrt_pred, harmful_fmt)
#    ws.write(row, col_idx_mg.IDX_MT, content_rec.mt, dflt_cell_fmt)
#    if content_rec.mt_harmful:
#        harmful_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_HARMFUL]]
#    else:
#        harmful_fmt = dflt_cell_fmt
#    ws.write(row, col_idx_mg.IDX_MTPRED, content_rec.mt_pred, harmful_fmt)
#    # write 'others' information
#    others_col_idx = col_idx_mg.IDX_MTPRED
#    for item in content_rec.others:
#        others_col_idx += 1
#        ws.write(row, others_col_idx, item, dflt_cell_fmt)
#    # get cell format and write zygosities
#    zygo_col_idx = content_rec.n_master_cols - 1
#    for pat_zygo in content_rec.pat_zygos:
#        zygo_col_idx += 1
#        if rare and content_rec.all_mutated:
#            zygo_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_HARMFUL]]
#        elif pat_zygo.shared_mutation and pat_zygo.is_hom:
#            zygo_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_HOM_SHARED]]
#        elif pat_zygo.shared_mutation:
#            zygo_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_SHARED]]
#        elif rare and pat_zygo.is_mutated:
#            zygo_fmt = cell_fmt_mg.cell_fmts[cell_colors[CELL_TYPE_RARE]]
##        elif pat_zygo.is_mutated:
##            zygo_fmt = cell_fmt
#        else:
#            zygo_fmt = dflt_cell_fmt
#        ws.write(row, zygo_col_idx, pat_zygo.zygo, zygo_fmt)
#    # add extra attributes
#    for attrib_idx in xrange(len(xtra_attribs)):
#        attrib = xtra_attribs[attrib_idx]
#        cell_attrib = None
#        if attrib == ATTRIB_RARE:
#            if rare:
#                cell_attrib = "yes"
#            else:
#                cell_attrib = "no"
#        if attrib == ATTRIB_HAS_SHARED:
#            if content_rec.has_shared_mutation:
#                cell_attrib = "yes"
#            else:
#                cell_attrib = "no"
#        if attrib == ATTRIB_STUDY:
#            if marked_color is not None:
#                cell_attrib = "yes"
#            else:
#                cell_attrib = "no"
#        if attrib == ATTRIB_CASES_GE_CTRLS:
#            if content_rec.cases_ge_ctrls:
#                cell_attrib = "yes"
#            else:
#                cell_attrib = "no"
#        if attrib == ATTRIB_HAS_MUTATION:
#            if content_rec.has_mutation:
#                cell_attrib = "yes"
#            else:
#                cell_attrib = "no"
#        ws.write(row, rec_size+attrib_idx, cell_attrib, dflt_cell_fmt)
#        
#def add_sheet(wb, sheet_name):
#    ws = wb.add_worksheet(sheet_name)
#    ws.set_default_row(12)
#    return ws
#
#def add_muts_sheet(wb, cell_fmt_mg, muts_rep, xtra_attribs):
#    ws = add_sheet(wb, muts_rep.sheet_name)
#    #ws = wb.add_worksheet(muts_rep.sheet_name)
#    ws.set_default_row(12)
#    mut_rec_size = muts_rep.record_size
#    write_header(ws,
#                 cell_fmt_mg,
#                 muts_rep.header_rec,
#                 mut_rec_size,
#                 muts_rep.col_idx_mg,
#                 xtra_attribs)
#    # write content
#    row = 1
#    for mut_rec in muts_rep.mut_recs:
#        including = True
#        for criterian in inc_criteria:
#            if (criterian == INC_SHARED_MUTATION) and ( not mut_rec.has_shared_mutation):
#                including = False
#                break
#        if including:
#            write_content(ws,
#                          cell_fmt_mg,
#                          row,
#                          mut_rec,
#                          mut_rec_size,
#                          muts_rep.col_idx_mg,
#                          xtra_attribs)
#            row += 1
#    set_layout(ws, mut_rec_size+len(xtra_attribs), muts_rep.col_idx_mg) 
#        
#def add_addn_csv_sheet(wb, dflt_cell_fmt, sheet_name, csv_file):
#    ws = add_sheet(wb, sheet_name)
#    with open(csv_file, 'rb') as csvfile:
#        csv_recs = list(csv.reader(csvfile, delimiter='\t'))
#        csv_row = 0
#        for xls_row in xrange(len(csv_recs)):
#            csv_rec = csv_recs[xls_row]
#            for col in xrange(len(csv_rec)):
#                ws.write(csv_row, col, csv_rec[col], dflt_cell_fmt)
#            csv_row += 1
#        csvfile.close()
#    ws.freeze_panes(1, 0)
#
## ****************************** main codes ******************************
#new_section_txt(" Generating report ")
#
#wb = xlsxwriter.Workbook(out_file)
#cell_fmt_mg = CellFormatManager(wb, COLOR_RGB)
#dflt_cell_fmt = cell_fmt_mg.default_format
#debug(cell_fmt_mg)
#
#for main_csv in csvs_list:
#    (sheet_name, sheet_csv) = main_csv.split(',')
#    muts_rep = MutationsReport(file_name=sheet_csv,
#                               n_master_cols=n_master_cols,
#                               sheet_name=sheet_name,
#                               color_region_infos=color_region_infos,
#                               freq_ratios=frequency_ratios)
#    debug(muts_rep)
#    mylogger.info("adding mutations sheet: " + sheet_name)
#    add_muts_sheet(wb, cell_fmt_mg, muts_rep, xtra_attribs)
#
#for addn_csv in addn_csvs_list:
#    (sheet_name, sheet_csv) = addn_csv.split(',')
#    add_addn_csv_sheet(wb, dflt_cell_fmt, sheet_name, sheet_csv)
#
#wb.close()
#
#new_section_txt(" F I N I S H <" + script_name + "> ")
