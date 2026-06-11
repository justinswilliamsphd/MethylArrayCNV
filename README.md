# MethylArrayCNV
Process and integrate DNAm array data for single-sample and cohort-level CNV analysis.
**NOTE**: The default references provided are of normal brain tissue and mouse brain tissue, however a custom reference set is easily implemented by providing those samples/GSM IDs in the sample sheet. Samples from GEO are downloaded locally.

Generate CNV profiles:
1) cohort level karyotype plots
2) sample-by-sample detailed plots
3) focal detailed annotation plots

optional:
1) create a normal reference set (saved as Rdata) 
2) create a generic annotation (focal analysis, saved as Rdata)
3) Include samples from GEO by GSM ID (downloads locally) 

## Environment requirements:
- r-base=4.4.3
- r-argparse=2.2.5
- r-plyr=1.8.9
- r-dplyr=1.1.4
- r-stringr=1.5.1
- bioconductor-minfi=1.52.0
- bioconductor-sesame=1.24.0
- bioconductor-sesamedata=1.24.0
- conumee2_2.1.2
- GEOquery_2.74.0
- bioconductor-rtracklayer=1.66.0
- GenomeInfoDb_1.42.0 
- ComplexHeatmap_2.22.0 
- r-circlize=0.4.16
- plyranges_1.22.0
- biomaRt_2.62.0
- IlluminaHumanMethylation450kmanifest_0.4.0
- IlluminaHumanMethylationEPICmanifest_0.3.0
- IlluminaHumanMethylationEPICv2manifest_1.0.0
- IlluminaHumanMethylationEPICv2anno.20a1.hg38_1.0.0
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
- IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1 
- IlluminaMouseMethylationanno.12.v1.mm10_0.0.2
- IlluminaMouseMethylationmanifest_0.0.1

## To Run:
### required parameters:  

&emsp;sample_sheet, full path to tab delimited sample name/ IDAT matchup  

### optional parameters:  

&emsp;**--id**, additional identifier for output (default: "")  

&emsp;**--analysis**, standard, karyotype, or both (default: both)  

&emsp;**--ref**, RData CNV reference obj (default: NULL)  
&emsp;&emsp;**NOTE**: if provided, will override sample sheet reference samples
  
&emsp;**--annot**, RData CNV annotation obj	(default: NULL)	 

&emsp;**--rn**, If specified, rename samples from Sentrix to SampleID (default: TRUE)  

&emsp;**--mdr**, If specified, create default reference set (default: FALSE)  
&emsp;&emsp;**NOTE**: if provided, will override sample sheet reference samples  
  
&emsp;**--mdr_human_sheet**, The default reference sample_sheet path (default: "default_human_brain_reference_450K.txt")  

&emsp;**--mdr_mouse_sheet**, The default reference sample_sheet path (default: "default_mouse_brain_reference_450K.txt")  

&emsp;**--mda**, If specified, create default annotation set (default: TRUE)  

&emsp;**--xy**, If specified, include XY chromosomes (default: FALSE)  

&emsp;**--kt**, Log2Ratio absolute threshold for inclusion in karyotype view (default: 0.1)  

## Input SampleSheet

| Index | SampleID | BaseName | ArrayType | GSMID | SampleType | Download | Annotation | AnnotColor |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | ---: |
| sentrix_id | A sample name | full_path_to_sentrix_id | version of array[450k, EPIC, EPICv2, mouse] | Accession from GEO | Query or Reference[query/ref]  | boolean[TRUE/FALSE] | Sample Grouping | Grouping Color |
| 205091240115_R05C01 | Sample1 | /full/path/205091240115_R05C01 | EPICv2 |  | query | FALSE | CancerType1 | #EDA909 |
|  | CONTR_ADENOPIT_02 |  | 450k | GSM2405408 | ref | TRUE |  |  |
|  | GBM_G34_sample1  |  | 450k | GSM2402854 | query | TRUE | GBM_G34 | #4287f5 |
