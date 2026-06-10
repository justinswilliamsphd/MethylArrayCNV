#------------------------------------------------------------------------------
# Generate CNV Profiles
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
get_line_args <- function(){
  check_dependencies(c("argparse"), FALSE)
	parser <- argparse::ArgumentParser(description='Create CNV profiles from methylation')
	parser$add_argument('sample_sheet', type="character", 
						help='Tab delimited sample name/ IDAT matchup')
  parser$add_argument('--id', default="", type="character", 
            help='Additional identifier for output (default: "")')
  parser$add_argument('--analysis', default="both", type="character", 
            help='standard, karyotype, or both (default: both)')
	parser$add_argument('--ref', default=NULL, type="character", 
						help='RData CNV reference obj')
	parser$add_argument('--annot', default=NULL, type="character", 
						help='RData CNV annotation obj')		
  parser$add_argument('--rn', default=TRUE, type="logical", 
            help='If specified, rename samples from Sentrix to SampleID')
	parser$add_argument('--mdr', default=FALSE, type="logical", 
						help='If specified, create default reference set')
  mdr_human_sheet=file.path(getwd(),"default_human_brain_reference_450K.txt")
  parser$add_argument('--mdr_human_sheet', default=mdr_human_sheet, type="character", 
            help='The default human reference sample_sheet path')
  mdr_mouse_sheet=file.path(getwd(),"default_mouse_brain_reference_Mouse.txt")
  parser$add_argument('--mdr_mouse_sheet', default=mdr_mouse_sheet, type="character", 
            help='The default mouse reference sample_sheet path')
	parser$add_argument('--mda', default=TRUE, type="logical", 
						help='If specified, create default annotation set (default: TRUE)')
  parser$add_argument('--xy', default=FALSE, type="logical", 
            help='If specified, include XY chromosomes (default: FALSE)')
  parser$add_argument('--kt', default=0.1, type="numeric", 
            help='Log2Ratio abs threshold for inclusion in karyotype (default: 0.1)')
	return(parser$parse_args())
}

#------------------------------------------------------------------------------
check_dependencies<-function(library_list=c(), load_packages=TRUE, load_manifest="ALL"){
  if((load_packages)&(length(library_list)==0)){
    message("Checking standard dependency requirements")    
  }else if((load_packages)&(length(library_list)>0)){
    message(paste("Loading standard packages, and checking dependencies:",paste0(library_list,collapse=","))) 
  }else{
    message("Checking/loading standard dependencies")    
  }
  #-------
  if(length(library_list)==0){
    library_list <- c(
      'argparse', 'plyr', 'dplyr', 'stringr', 'minfi', 'sesame', 'conumee2',
      'GEOquery', 'rtracklayer', 'GenomeInfoDb', 'ComplexHeatmap',
      'circlize', 'plyranges',
      'IlluminaHumanMethylation450kmanifest',
      'IlluminaHumanMethylationEPICmanifest', 
      'IlluminaHumanMethylationEPICv2manifest',
      'IlluminaHumanMethylationEPICv2anno.20a1.hg38',
      'IlluminaHumanMethylationEPICanno.ilm10b4.hg19',
      'IlluminaHumanMethylation450kanno.ilmn12.hg19',
      'IlluminaMouseMethylationanno.12.v1.mm10',
      'IlluminaMouseMethylationmanifest',
      "biomaRt"
    )    
  }
  test_list <- c()
  for(lib in library_list){
    lib_test <- try(find.package(lib), silent = TRUE)
    if("try-error" == class(lib_test)) {
      test_list <- c(lib, test_list)
    }
  }
  if(length(test_list) > 0) {
    test_str=paste(test_list, collapse=", ")
    warning = paste("Packages:", test_str, " are not available")
    stop(warning)
  }
  # light weight packages that are easier to load:
  if(load_packages){
    suppressPackageStartupMessages(library(plyr, quietly = TRUE))
    suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
    suppressPackageStartupMessages(library(conumee2, quietly = TRUE))
    if(load_manifest=="ALL"){
      suppressPackageStartupMessages(library(IlluminaMouseMethylationanno.12.v1.mm10, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaMouseMethylationmanifest, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaHumanMethylation450kmanifest, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICmanifest, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, quietly = TRUE))  
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38, quietly = TRUE))    
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2manifest, quietly = TRUE))      
    }
    if("mouse" %in% load_manifest){
      suppressPackageStartupMessages(library(IlluminaMouseMethylationanno.12.v1.mm10, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaMouseMethylationmanifest, quietly = TRUE))
    }
    if("450k" %in% load_manifest){
      suppressPackageStartupMessages(library(IlluminaHumanMethylation450kmanifest, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19, quietly = TRUE))
    }    
    if("EPIC" %in% load_manifest){
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICmanifest, quietly = TRUE))
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, quietly = TRUE))  
    }   
    if("EPICv2" %in% load_manifest){
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38, quietly = TRUE))    
      suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2manifest, quietly = TRUE))
    }   
  }
  message("-dependency requirements met")
  return(TRUE)
}

#------------------------------------------------------------------------------
log_message <- function(log_file, msg=msg, append=TRUE){
  msg = paste0("\n", msg)
  cat(msg)
  cat(file=log_file, text=msg, append=append)
}

#------------------------------------------------------------------------------
check_file_exists<-function(line_args, file_list){
  for (filename in names(file_list)){
    file_path = file_list[[filename]]
    if(!is.null(file_path)){
      if(!file.exists(file_path)){
        warning = paste0(filename,": ",file_path,
          "\ndoes not exist. Exiting...")
        log_message(line_args$log_file, msg=warning)
        stop(warning)
      }else{
        msg = paste0(filename,": ",file_path)
        log_message(line_args$log_file, msg=msg)        
      }
    }else{
      msg = paste("Not using:",filename)
      log_message(line_args$log_file, msg=msg)   
    }
  }
}

#------------------------------------------------------------------------------
check_line_args <- function(line_args){
  # input paths:
  line_args$sample_sheet = normalizePath(line_args$sample_sheet)
  # output paths:
  timestamp = strtrim(gsub("[: ]", ".", Sys.time()), 16)
  if(line_args$id != ""){
    line_args$outdir = file.path(getwd(),paste(timestamp,line_args$id,"CNV_Analysis",sep="_"))
  }else{
    line_args$outdir = file.path(getwd(),paste0(timestamp,"_CNV_Analysis"))
  }
  dir.create(line_args$outdir, showWarnings=FALSE, recursive=TRUE)
  line_args$output_analysis_prefix = file.path(line_args$outdir, timestamp)
  # sample genome directory
  line_args$genome_plot_outdir=paste0(line_args$output_analysis_prefix,"_sample_genome_output")
  dir.create(line_args$genome_plot_outdir, showWarnings=FALSE, recursive=TRUE)
  # focal sample plots
  line_args$focal_output_outdir=paste0(line_args$output_analysis_prefix,"_sample_focal_output")
  dir.create(line_args$focal_output_outdir, showWarnings=FALSE, recursive=TRUE)
  # cohort directory
  line_args$cohort_output_outdir=paste0(line_args$output_analysis_prefix,"_cohort_genome_output")
  dir.create(line_args$cohort_output_outdir, showWarnings=FALSE, recursive=TRUE)
  # RData directory
  line_args$RData_outdir=paste0(line_args$output_analysis_prefix,"_RData")
  dir.create(line_args$RData_outdir, showWarnings=FALSE, recursive=TRUE)
  # create log file path
  line_args$log_file = paste0(
    line_args$output_analysis_prefix,
    "_run_log.txt"
  )
  msg = paste("Run started:", date())
  log_message(line_args$log_file, msg=msg, append=F)  
  # ensure samplesheet exists:
  check_files = list(
    sample_sheet = line_args$sample_sheet,
    ref = line_args$ref,
    annot = line_args$annot
  )
  check_file_exists(line_args, check_files)
  return(line_args)
}

#------------------------------------------------------------------------------
load_sample_sheet<-function(line_args){
  msg = paste0("sample_sheet: ", line_args$sample_sheet)
  log_message(line_args$log_file, msg=msg)
  ##  
	samples = read.delim(line_args$sample_sheet, stringsAsFactors = FALSE)

  missing_index = samples[is.na(samples$Index),]
  with_index = samples[is.na(samples$Index)==FALSE,]
  if(dim(missing_index)[1]>0){
    missing_GSM = missing_index[is.na(missing_index$GSMID),]
    if(dim(missing_GSM)[1]==0){
      missing_index$Index = missing_index$GSMID
      msg = paste0("missing indices: ", dim(missing_index)[1], " replaced with GSM")
      log_message(line_args$log_file, msg=msg)
      samples=rbind(with_index,missing_index)
    }else{
      msg = paste0("Samples with Index and GSM missing, please check sample sheet. Exiting...")
      log_message(line_args$log_file, msg=msg)
      stop()
    }
  }

	samples$Index=as.character(samples$Index)
	rownames(samples) = samples$Index

  samples$Download = as.logical(samples$Download) 

	print(head(samples))
	line_args$sample_sheet = samples
	return(line_args)
}

#------------------------------------------------------------------------------
rename_samples<-function(rgset, sample_sheet){ 
  # load and reorder sample names by rgset colname order
  rownames(sample_sheet) = sample_sheet$Index
  sample_sheet = sample_sheet[colnames(rgset),]
  colnames(rgset)<-sample_sheet$SampleID
  # deprecated: now we reference the samples by SampleID, so change the table index:
  #rownames(sample_sheet) = sample_sheet$SampleID
  #return_list = list(rgset=rgset, sample_sheet=sample_sheet)
  return(rgset) 
}

#------------------------------------------------------------------------------
get_gsm_basename<-function(idat_list){
  bn = stringr::str_replace(idat_list, "_Grn.idat", "")
  bn = unique(stringr::str_replace(bn, "_Red.idat", ""))
  return(bn)
}

#------------------------------------------------------------------------------
download_samples<-function(line_args, sample_sheet, dl_dir){
  for(gsm_id in sample_sheet$GSMID){
    dl_bool = sample_sheet[sample_sheet$GSMID==gsm_id,]$Download
    if(dl_bool==TRUE){
      GEOquery::getGEOSuppFiles(GEO=gsm_id, baseDir=dl_dir, makeDirectory = FALSE)
      gz_idats = list.files(path = dl_dir, pattern = paste0(gsm_id,".*.idat.gz"),
        full.names=TRUE)
      for(gz_idat in gz_idats){
        R.utils::gunzip(gz_idat)      
      }
      idat_list = list.files(path = dl_dir, pattern = paste0(gsm_id,".*.idat"),
        full.names=TRUE)
      print(idat_list)
      bn = get_gsm_basename(idat_list)
      index_name = stringr::str_split(bn,"/")
      index_name = index_name[[1]][length(index_name[[1]])]
      print(paste(index_name,":",bn))
      sample_sheet[sample_sheet$GSMID==gsm_id,"BaseName"] = bn
      sample_sheet[sample_sheet$GSMID==gsm_id,"Index"] = index_name
      msg = paste0("unpacked sample: ", gsm_id, ", ", bn)
      log_message(line_args$log_file, msg=msg)      
    }      
  }
  rownames(sample_sheet) = sample_sheet$Index
  return(sample_sheet)
}

#------------------------------------------------------------------------------
load_sample_data<-function(line_args, sample_sheet, sample_type){
  type_sample_sheet = sample_sheet[sample_sheet$SampleType==sample_type,] 
  array_mset_list = list()
  for(array_type in unique(type_sample_sheet$ArrayType)){
    msg = paste0("loading samples: ", sample_type," / ", array_type)
    log_message(line_args$log_file, msg=msg)
    ## 
    array_sample_sheet = type_sample_sheet[type_sample_sheet$ArrayType==array_type,]
    dl_array_sample_sheet = array_sample_sheet[array_sample_sheet$Download==TRUE,]
    if(dim(dl_array_sample_sheet)[1]>0){
      if(dim(dl_array_sample_sheet)[1] == dim(array_sample_sheet)[1]){
        msg = paste0("all samples require download: ", sample_type," / ", array_type)
        log_message(line_args$log_file, msg=msg)
        ##
        dl_dir = paste(line_args$output_analysis_prefix,sample_type,array_type,"downloads",sep="_")
        dir.create(dl_dir, showWarnings=FALSE, recursive=TRUE)
        array_sample_sheet = download_samples(line_args, array_sample_sheet, dl_dir)       
      }else if(dim(dl_array_sample_sheet)[1] < dim(array_sample_sheet)[1]){
        msg = paste0("partial dataset download: ", sample_type," / ", array_type)
        log_message(line_args$log_file, msg=msg)
        ##
        dl_dir = dirname(dl_array_sample_sheet$BaseName)
        array_sample_sheet = download_samples(line_args, array_sample_sheet, dl_dir)  
      }
    }
    # saving dataset sample_sheet to file
    msg = paste0("saving sample_sheet to file: ", sample_type,", ", array_type)
    log_message(line_args$log_file, msg=msg)
    sample_sheet_outpath = paste(
      line_args$output_analysis_prefix, sample_type, array_type, "sample_sheet.txt",  sep="_")
    write.table(array_sample_sheet, sample_sheet_outpath, sep = "\t")
    ##  
    msg = paste0("building rgset: ", sample_type,", ", array_type)
    log_message(line_args$log_file, msg=msg)
    # build RGset class
    rgset = minfi::read.metharray(array_sample_sheet$BaseName, force= TRUE, verbose = TRUE)
    if(array_type=="EPICv2"){
      annotation(rgset)["array"] = "IlluminaHumanMethylationEPICv2"
      annotation(rgset)["annotation"] = "20a1.hg38"      
    }
    if(array_type=="mouse"){
      annotation(rgset)["array"] = "IlluminaMouseMethylation"
      annotation(rgset)["annotation"] = "12.v1.mm10"      
    }
    ##  
    msg = paste0("checking sample set: ", sample_type,", ", array_type)
    log_message(line_args$log_file, msg=msg)
    # compare sample names from index and dataset:
    diff_samples = setdiff(colnames(rgset),rownames(array_sample_sheet))
    if(length(diff_samples) >= 1)
      {
        msg = paste0(
          "Provided Sample Index does not match dataset:\n",
          paste0(diff_samples,collapse=", "), "\n",
          "This could be benign, due to downloaded sample Index.\n",
          "If not, consider correcting and re-running."
        )
        log_message(line_args$log_file, msg=msg)
      }
    if((line_args$rn)&(sample_type != "ref")){
      msg = paste0("renaming samples: ", sample_type,", ", array_type)
      log_message(line_args$log_file, msg=msg)      
      # renaming samples
      rgset = rename_samples(rgset, array_sample_sheet)
    }
    mset = minfi::preprocessNoob(rgset, dyeMethod="single")
    array_mset_list[[array_type]] = mset
  }
  return(array_mset_list)
}

#------------------------------------------------------------------------------
infer_array_type<-function(mset){
  array_type = minfi::annotation(mset)[[1]]
  array_type = stringr::str_replace(array_type, "IlluminaHumanMethylation", "")
  return(array_type)
}

#------------------------------------------------------------------------------
grange_hg19_to_hg38<-function(line_args, hg19_gr){
  gz_dest_path = file.path(line_args$outdir,"hg19ToHg38.over.chain.gz")
  dest_path = file.path(line_args$outdir,"hg19ToHg38.over.chain")
  if(!file.exists(dest_path)){
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    download.file(url, gz_dest_path, method="curl")
    R.utils::gunzip(gz_dest_path, overwrite=TRUE) 
  }
  ch = rtracklayer::import.chain(dest_path)
  GenomeInfoDb::seqlevelsStyle(hg19_gr) = "UCSC"
  hg38_gr = unlist(as(rtracklayer::liftOver(hg19_gr, ch), "GRangesList"))  
  return(hg38_gr)  
}

#------------------------------------------------------------------------------
hs_mm_gene_coords<-function(gene_list, genome){
  # note genome is the genome of the gene symbols PROVIDED
  # if 'hs' is supplied, mouse geneSymbols will be returned
  #
  # load match-ups
  urlpath=paste0(
    "https://www.informatics.jax.org/downloads/reports/",
    "HOM_MouseHumanSequence.rpt"
  )
  gene_matchup= read.csv(urlpath, sep="\t")
  # split by genome
  hsg = gene_matchup[gene_matchup$Common.Organism.Name=="human",]
  mmg = gene_matchup[gene_matchup$Common.Organism.Name=="mouse, laboratory",]
  converted_gene_list = c()
  if(genome=="hs"){
    rownames(hsg) = make.names(hsg$Symbol, unique=TRUE)
    rownames(mmg) = make.names(mmg$DB.Class.Key,unique=TRUE)
    found = intersect(gene_list,  rownames(hsg))
    idx = make.names(hsg[found,'DB.Class.Key'])
    converted_gene_list = mmg[idx,"Symbol"]
    mart = biomaRt::useMart(
      "ensembl", 
      dataset = "mmusculus_gene_ensembl"
    )
  } else if (genome=="mm"){
    rownames(mmg) = make.names(mmg$Symbol, unique=TRUE)
    rownames(hsg) = make.names(hsg$DB.Class.Key,unique=TRUE)
    found = intersect(gene_list,  rownames(mmg))
    idx = make.names(mmg[found,'DB.Class.Key'])
    converted_gene_list = hsg[idx,"Symbol"]
    mart = biomaRt::useMart(
      "ensembl", 
      dataset = "hsapiens_gene_ensembl"
    )    
  } else {
    warning = paste("Genome must be hs or mm. Exiting")
    #log_message(line_args$log_file, msg=warning)
    stop(warning)  
  }
  print(paste0("recovered genes: ",length(converted_gene_list), "/", length(gene_list)))

  g_df = biomaRt::getBM(
    attributes = c("external_gene_name", 
                    "chromosome_name", 
                    "start_position", 
                    "end_position"),
    filters = "external_gene_name",
    values = converted_gene_list,
    mart = mart
  )
  colnames(g_df) = c("name","chr","start","end")
  g_df$chr = paste0("chr",g_df$chr) 
  gr = GenomicRanges::makeGRangesFromDataFrame(
      g_df, keep.extra.columns=TRUE
  )
  names(gr) = gr$name
  return(gr)
}

#------------------------------------------------------------------------------
make_detail_gr<-function(detail_regions, cancer_genes){
  # format cancer genes
  cancer_genes$name = cancer_genes$SYMBOL
  mcols(cancer_genes)[c("SYMBOL")]=NULL
  # format detail_regions
  names(detail_regions) = detail_regions$name
  del_cols = c("score", "thick", "probes_gene", "probes_promoter")
  del_cols = intersect(names(mcols(detail_regions)),del_cols)
  mcols(detail_regions)[del_cols]=NULL
  # intersect and combine datasets
  # build a cancer-gene compatible detail regions
  common_genes = intersect(cancer_genes$name, detail_regions$name)
  detail_only = setdiff(detail_regions$name,cancer_genes$name)
  select_cc_only = intersect(c("SMARCB1","Smarcb1"),names(cancer_genes))
  glist=list()
  glist[[1]] = detail_regions[detail_only]
  glist[[2]]  = cancer_genes[common_genes]
  glist[[3]]  = cancer_genes[select_cc_only]
  update_detail_regions = do.call(c, as(glist, "GRangesList"))
  return(update_detail_regions)
}

#------------------------------------------------------------------------------
make_annots<-function(line_args, all_array_types, karyo_annot_set=NULL){
  msg = paste0("creating annotations for arrays: ", paste0(all_array_types, collapse=", "))
  log_message(line_args$log_file, msg=msg)
  cat("")
  # conditionals: EPICv2 is the only hg38
  # if any hg19 are used, defaults to hg19 (by design: conumee2)

  if(!is.null(karyo_annot_set)&(all("mouse" != all_array_types))){
    genome = karyo_annot_set@genome
    genome$strand = NA
    # p
    genome_p = genome
    genome_p$name = paste0(genome_p$chr,"p")
    genome_p$start = 1
    genome_p$end = genome_p$pq-1
    # q
    genome_q = genome
    genome_q$name = paste0(genome_q$chr,"q")
    genome_q$start = genome_q$pq
    genome_q$end = genome_q$size
    # combine
    g_df = rbind(genome_p,genome_q)
    g_df = g_df[,c("chr","start","end","strand","name")]
    rownames(g_df) = g_df$name
    detail_regions = GenomicRanges::makeGRangesFromDataFrame(
      g_df, keep.extra.columns=TRUE
    )
  }else if(!is.null(karyo_annot_set)&(all("mouse" == all_array_types))){
    genome = karyo_annot_set@genome
    colnames(genome) = c("chr", "end")
    genome$strand = NA
    genome$start = 1
    genome$name = rownames(genome)
    genome = genome[,c("chr","start","end","strand","name")]
    detail_regions = GenomicRanges::makeGRangesFromDataFrame(
      genome, keep.extra.columns=TRUE
    )
  }else if(all("mouse" == all_array_types)){
    # NOTE: mouse arays not eligible for focal analysis.
    #
    data(detail_regions.mm10)
    mcols(detail_regions)[c("thick")] = NULL
    # create an enhanced detail_regions
    cancer_genes_list = c("PPM1D", "TERT", "OTX2", "GFI1", "GFI1B", "CCND2", "GLI1", "PRDM6")
    cancer_gene_gr = hs_mm_gene_coords(cancer_genes_list, "hs")
    detail_regions = c(detail_regions, cancer_gene_gr)
    
  }else if(all("EPICv2" == all_array_types)==FALSE){
    # hg19
    data(exclude_regions)
    data(detail_regions)
    data(cancer_genes_hg19)
    # add OTX2 hg19 chr14:57,267,425-57,277,197
    OTX2 = GenomicRanges::GRanges("chr14",IRanges(57267425,57277197),SYMBOL="OTX2")
    names(OTX2) = c("OTX2")
    cancer_genes = c(cancer_genes, OTX2)
    detail_regions = make_detail_gr(detail_regions, cancer_genes)
    # cancer_genes formatting
    cancer_genes$name = cancer_genes$SYMBOL
  }else{
    # hg38
    data(exclude_regions)
    data(detail_regions.hg38)
    data(cancer_genes_hg38)
    # add OTX2 hg38 chr14:56,799,905-56,816,693
    OTX2 = GenomicRanges::GRanges("chr14",IRanges(56799905,56816693),SYMBOL="OTX2")
    names(OTX2) = c("OTX2")
    cancer_genes = c(cancer_genes, OTX2)
    detail_regions = make_detail_gr(detail_regions, cancer_genes)
    exclude_regions = grange_hg19_to_hg38(line_args, exclude_regions)
    #
    cancer_genes$name = cancer_genes$SYMBOL
  }
  if(!line_args$xy!=TRUE){
    if(all("mouse" == all_array_types)){
      exclude_regions=NULL
      detail_regions = detail_regions[!(seqnames(detail_regions) %in% c("chrX", "chrY"))]
      cancer_genes = NULL
    }else{
      exclude_regions = exclude_regions[!(seqnames(exclude_regions) %in% c("chrX", "chrY"))]
      detail_regions = detail_regions[!(seqnames(detail_regions) %in% c("chrX", "chrY"))]
      cancer_genes = cancer_genes[!(seqnames(cancer_genes) %in% c("chrX", "chrY"))]      
    }
  }
  anno = conumee2::CNV.create_anno(
    array_type = all_array_types, 
    exclude_regions = if(all("mouse" == all_array_types)){NULL}else{exclude_regions}, 
    detail_regions = detail_regions,
    chrXY=line_args$xy)
  return(anno)
}

#------------------------------------------------------------------------------
core_conumee_analysis<-function(line_args, all_array_types, ref_mset_list, query_mset_list, annot_set){
  msg = "normalizing and combining intensity values"
  log_message(line_args$log_file, msg=msg)     
  # combine intensity values
  ref_cnv_obj = conumee2::CNV.load(ref_mset_list[[1]])
  fit_cnv_list=c()
  for(query_array_type_num in 1:length(query_mset_list)){
    query_cnv_obj = conumee2::CNV.load(query_mset_list[[query_array_type_num]])
    fit_cnv_i=conumee2::CNV.fit(query_cnv_obj, ref_cnv_obj, annot_set)
    if(query_array_type_num==1){
      fit_cnv=fit_cnv_i
    }else if(query_array_type_num>1){
      fit_cnv = conumee2::CNV.combine(fit_cnv_i, fit_cnv)
    }
  }
  msg = "binning probes (provided by annotations)"
  log_message(line_args$log_file, msg=msg)    
  fit_cnv = conumee2::CNV.bin(fit_cnv)
  
  msg = "Adding gene definitions"
  log_message(line_args$log_file, msg=msg)    
  fit_cnv = conumee2::CNV.detail(fit_cnv)
  
  msg = "Segment genome by copy-number state"
  log_message(line_args$log_file, msg=msg)    
  fit_cnv = conumee2::CNV.segment(fit_cnv)

  if(all("mouse" == all_array_types)){
    msg = "Focal analysis cannot be run on mouse arrays"
    log_message(line_args$log_file, msg=msg)
  }else{
    msg = "Performing focal CNV analysis"
    log_message(line_args$log_file, msg=msg)
    consensus_cancer_genes_hg19=cancer_genes
    consensus_cancer_genes_hg38=cancer_genes 
    fit_cnv = conumee2::CNV.focal(
      fit_cnv,
      sig_cgenes=TRUE
    )  
  }
  return(fit_cnv)
} 

#------------------------------------------------------------------------------
concat_func<-function(input){
  not_blank_input = input[input!=""] 
  unique_input = unique(not_blank_input)
  paste(unique_input, collapse=";")
}

#------------------------------------------------------------------------------
focal_report<-function(line_args, fit_cnv){
  # make a focal report
  focal = conumee2::CNV.write(fit_cnv, what="focal")

  focal_genes = list()
  focal_genes$amp_details = as.data.frame(dplyr::bind_rows(focal$"amplified detail regions"))
  rownames(focal_genes$amp_details) = names(focal$"amplified detail regions")
  
  focal_genes$del_details = as.data.frame(dplyr::bind_rows(focal$"deleted detail regions"))
  rownames(focal_genes$del_details) = names(focal$"deleted detail regions")

  focal_genes$amp_cc_genes = as.data.frame(dplyr::bind_rows(focal$"amplified genes from the Cancer Gene Census"))
  rownames(focal_genes$amp_cc_genes) = names(focal$"amplified genes from the Cancer Gene Census")

  focal_genes$del_cc_genes = as.data.frame(dplyr::bind_rows(focal$"deleted genes from the Cancer Gene Census"))
  rownames(focal_genes$del_cc_genes) = names(focal$"deleted genes from the Cancer Gene Census")

  combined = t(dplyr::bind_cols(focal_genes,.name_repair="minimal")) 
  combined_df = as.data.frame(combined, check.names=FALSE, make.names=NA)
  combined_df$Genes = rownames(combined)
  combined_df[is.na(combined_df)]=""

  focal_report = plyr::ddply(combined_df, .(Genes), colwise(concat_func))
  focal_report = as.data.frame(focal_report)
  rownames(focal_report) = focal_report$Genes
  focal_report$Genes = NULL
  outfile = file.path(line_args$focal_output_outdir,
      "standard_analysis_cohort_FOCAL_CNV.txt")
  write.table(focal_report, outfile, sep="\t", quote=FALSE)
}

#------------------------------------------------------------------------------
standard_conumee_analysis<-function(line_args, all_array_types, fit_cnv){
  msg = "Standard conumee analysis"
  log_message(line_args$log_file, msg=msg)  
  #  
  if(all("mouse" == all_array_types)){
    msg = "Skipping focal analysis report"
    log_message(line_args$log_file, msg=msg)
  }else{
    msg = "Creating focal analysis report"
    log_message(line_args$log_file, msg=msg)
    focal_report(line_args, fit_cnv)
  }

  msg = "Plotting genome CNV profiles by sample\n"
  log_message(line_args$log_file, msg=msg)    
  conumee2::CNV.genomeplot(
    fit_cnv,
    detail=TRUE,
    sig_cgenes=TRUE,
    nsig_cgenes=3,
    output="pdf",
    directory=line_args$genome_plot_outdir
  )

  if(all("mouse" == all_array_types)){
    msg = "Skipping detailed regions report"
    log_message(line_args$log_file, msg=msg)
  }else{
    msg = "Plotting all detailed regions"
    log_message(line_args$log_file, msg=msg)    
    conumee2::CNV.detailplot_wrap(
      fit_cnv,
      output="pdf",
      directory=line_args$focal_output_outdir
    )
  }
  msg = "Plotting summary level plots"
  log_message(line_args$log_file, msg=msg) 
  conumee2::CNV.summaryplot(
    fit_cnv,
    output="pdf",
    directory=line_args$cohort_output_outdir
  )
  # calculate height by number of samples:
  height = 1+(length(fit_cnv@fit$noise) * .25)
  conumee2::CNV.heatmap(
    fit_cnv,
    set_par=TRUE,
    output="pdf",
    directory=line_args$cohort_output_outdir,
    width=6,
    height=height
  )

  output_list = c("detail", "bins")
  file_exts = list(detail=".txt", bins=".igv", segments=".seg")
  for(output in output_list){
    outfile = file.path(line_args$cohort_output_outdir,
      paste0("standard_analysis_cohort_", output, file_exts[output]))
    conumee2::CNV.write(fit_cnv, file=outfile, what=output)    
  }
}

#------------------------------------------------------------------------------
karyotype_conumee_analysis<-function(line_args, all_array_types, fit_cnv){
  msg = "Saving karyotype CNV object"
  log_message(line_args$log_file, msg=msg)
  # format detailed region DF:
  details_list = conumee2::CNV.write(fit_cnv, what="detail")
  details_df = do.call(rbind, details_list)[,c("Name", "Sample", "Value")]
  details_df2 = reshape(details_df, direction="wide", idvar="Name", timevar="Sample")
  colnames(details_df2) = stringr::str_replace(colnames(details_df2), "Value.", "")
  details_df2 = details_df2[!is.na(details_df2$Name),]
  rownames(details_df2) = details_df2$Name
  details_df2$Name = NULL
  # for saving to file, include NA vals:
  if((line_args$xy==TRUE)&(all(all_array_types!="mouse"))){
    chroms = paste0("chr",c(1:22,"X","Y"))
  }else if((line_args$xy==TRUE)&(all(all_array_types=="mouse"))){
    chroms = paste0("chr",c(1:19,"X","Y"))
  }else if((line_args$xy==FALSE)&(all(all_array_types!="mouse"))){
    chroms = paste0("chr", 1:22)
  } else if((line_args$xy==FALSE)&(all(all_array_types=="mouse"))){
    chroms = paste0("chr", 1:19)
  }
  if(all(all_array_types!="mouse")){
    chr_arms =  plyranges::interweave(paste0(chroms,"p"), paste0(chroms,"q"))
    missing = setdiff(chr_arms,rownames(details_df2))
    details_df2[missing,]=NA
    details_df2 = details_df2[chr_arms,]
  }else{
    missing = c()
  }
  # save full
  outfile = file.path(line_args$cohort_output_outdir,
    "full_karyotype_cohort_details.txt")
  write.table(details_df2, file=outfile, sep="\t", quote=FALSE, 
    col.names=TRUE, row.names=FALSE)
  # remove NA for filtered output
  details_df2 = details_df2[!(rownames(details_df2) %in% missing),]
  # apply max value to show chr in cohort
  #l2r_idx = abs(apply(details_df2, 1, FUN = max))>line_args$kt
  #details_df2 = details_df2[l2r_idx,]
  # save to file
  #outfile = file.path(line_args$cohort_output_outdir,
  #  "filtered_karyotype_cohort_details.txt")
  #write.table(details_df2, file=outfile, sep="\t", quote=FALSE, 
  #  col.names=TRUE, row.names=FALSE)
  make_heatmap(line_args, details_df2)
}

#------------------------------------------------------------------------------
make_heatmap<-function(line_args, df){
  # parse the samplesheet for annotations
  sample_sheet = line_args$sample_sheet[line_args$sample_sheet$SampleType=="query",]
  if(line_args$rn==TRUE){
    rownames(sample_sheet) = sample_sheet$SampleID 
  }else{
    rownames(sample_sheet) = sample_sheet$Index    
  }
  # annotation color list
  color_list = setNames(as.character(sample_sheet$AnnotColor), as.character(sample_sheet$Annotations))
  color_list = color_list[!duplicated(names(color_list))]
  # create bottom column annotations
  plot_annots = ComplexHeatmap::HeatmapAnnotation(
    which="column",
    SampleType=sample_sheet[colnames(df),]$Annotations,
    col=list(SampleType=color_list),
    show_legend=TRUE,
    annotation_name_gp=gpar(fontsize = 8),
    height=unit(.125, "inches"),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 8),
      labels_gp=gpar(fontsize = 8))
  )
  # heatmap colors:
  legend_colors = c("darkblue", "white", "#F16729")
  legend_range = c(-0.25, 0, 0.25)
  legend_markers= c(-0.25, 0, 0.25)
  legend_labels = c("-0.25", "0", "0.25")
  mat_color_func = circlize::colorRamp2(legend_range, legend_colors)
  legend_matrix = ComplexHeatmap::Legend(
    title = "Log2Ratio", 
    col_fun = mat_color_func, 
    at = legend_markers, 
    labels = legend_labels
  )  
  #----

  # determine plot dimensions, geared toward fontsize 8
  body_width=(dim(df)[2]*.125)
  body_height=(dim(df)[1]*.125)
  fig_width = body_width +(max(nchar(rownames(df))+10)*.05) + 3
  fig_height = body_height + (max(nchar(colnames(df)))*.05) + 2

  karyotype_plotfile = file.path(line_args$cohort_output_outdir,"karyotype_heatmap.pdf")
  pdf(karyotype_plotfile, height=fig_height, width=fig_width)  
    hm = ComplexHeatmap::Heatmap(
      df,
      na_col = "white",
      col=mat_color_func,
      # dendrogram args
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      column_dend_reorder=TRUE,
      # sample naming
      row_title=NULL,
      column_title=NULL,
      show_row_names = TRUE,
      row_names_side = "left",
      show_column_names = TRUE,
      column_names_side = "bottom",
      #row_names_centered=TRUE,
      row_names_gp = gpar(fontsize = 8, just="left"),# was 10
      column_names_gp = gpar(fontsize = 8),# was 10
      column_names_rot = 90,
      # row, col formatting
      row_split=1:dim(df)[1],
      column_split=1:dim(df)[2],
      gap = unit(.025, "inches"),
      # annotation args
      bottom_annotation=plot_annots,
      # heatmap legends
      show_heatmap_legend = TRUE,
      heatmap_legend_param = list(
        title = "Log2Ratio", 
        title_gp=gpar(fontsize = 8),
        labels_gp=gpar(fontsize = 8),
        legend_gp=gpar(fontsize = 8)),
      # heatmap sizing
      #heatmap_width=unit(fig_width, "inches"),# whole heatmap
      width=unit(body_width, "inches"),# heatmap body
      #heatmap_height=unit(fig_height, "inches"),# whole heatmap
      height=unit(body_height, "inches")# heatmap body
    )
    ComplexHeatmap::draw(hm, padding = unit(c(1, 1, 1, 1), "inches"))
  dev.off()
}

#------------------------------------------------------------------------------
main<-function(){
  # return command line arguments
  line_args = get_line_args()
	# check for library availability:
  load_manifest = unique(read.delim(line_args$sample_sheet)$ArrayType)
	check_dependencies(load_manifest=load_manifest)
  #line_args = make_test_line_args()
	# check input, add output directories/prefixes to line_args
	line_args = check_line_args(line_args)
	# format the input sample list
	line_args = load_sample_sheet(line_args)
  line_args_path = file.path(line_args$RData_outdir,"line_args.RData")
  save(line_args, file=line_args_path)
  # load query data from samplesheet
  query_mset_list = load_sample_data(line_args, line_args$sample_sheet, "query")
  def_query_path = file.path(line_args$RData_outdir,"query_set.RData")
  save(query_mset_list, file=def_query_path)
  # load reference mset:
  if(!is.null(line_args$ref)){
    ref_mset_list = load(line_args$ref_set)
  }else if(line_args$mdr==TRUE){
    if("mouse" %in% line_args$sample_sheet$ArrayType){
      mdr_sheet = read.delim(line_args$mdr_mouse_sheet)
    }else{
      mdr_sheet = read.delim(line_args$mdr_human_sheet)
    }
    ref_mset_list = load_sample_data(line_args, mdr_sheet, "ref")
    def_ref_path = file.path(line_args$RData_outdir,
      paste0("default_reference_set.", paste(names(ref_mset_list),collapse="_"), ".RData"))
    save(ref_mset_list, file=def_ref_path)
  }else{
    ref_mset_list = load_sample_data(line_args, line_args$sample_sheet, "ref")
    cust_ref_path = file.path(line_args$RData_outdir,
      paste0("custom_reference_set.", paste(names(ref_mset_list),collapse="_"), ".RData"))
    save(ref_mset_list, file=cust_ref_path)
  }
  # create or load annotations 
  if(line_args$mda){
    ref_array_type = names(ref_mset_list)
    query_array_types = names(query_mset_list)
    all_array_types = unique(c(query_array_types, ref_array_type))
    annot_set = make_annots(line_args, all_array_types)
    # save to file
    def_annot_path = file.path(line_args$RData_outdir,
      "default_annotations.RData")
    save(annot_set, file=def_annot_path)    
  }else{
    annot_set = load(line_args$annot) 
  }
  msg = paste0("Starting analysis: ", line_args$analysis)
  log_message(line_args$log_file, msg=msg) 
  if(line_args$analysis %in% c("standard", "both")){
    fit_cnv = core_conumee_analysis(line_args, all_array_types, ref_mset_list, query_mset_list, annot_set)
    # save obj to file
    fit_cnv_path = file.path(line_args$RData_outdir,
        "standard_analysis_fit_cnv_obj.RData")
    save(fit_cnv, file = fit_cnv_path)
    standard_conumee_analysis(line_args, all_array_types, fit_cnv)
  }
  if(line_args$analysis %in% c("karyotype", "both")){
    msg = "creating karyotype annotation set"
    log_message(line_args$log_file, msg=msg)     
    karyo_annot_set = make_annots(line_args, all_array_types, annot_set)
    fit_cnv = core_conumee_analysis(line_args, all_array_types, ref_mset_list, query_mset_list, karyo_annot_set)
    # save this object
    fit_cnv_path = file.path(line_args$RData_outdir,
        "karyotype_analysis_fit_cnv_obj.RData")
    save(fit_cnv, file = fit_cnv_path)    
    karyotype_conumee_analysis(line_args, all_array_types, fit_cnv)
  }
  msg = "CNV analysis complete"
  log_message(line_args$log_file, msg=msg) 
}

#------------------------------------------------------------------------------
if (!interactive()) {
  #options(warn=-1)
  main()
  #options(warn=0)
}