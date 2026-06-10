# MethylArrayCNV
Process and integrate DNAm array data for single-sample and cohort-level CNV analysis.

Generate CNV profiles:
  1) cohort level karyotype plots
  2) sample-by-sample detailed plots
  3) focal detailed annotation plots

optional:
  4) create a normal reference Rdata for CNS/Brain samples
  5) create a generic annotation (focal) RData file

## Environment requirements:
1) r-base=4.4.3
2) r-argparse=2.2.5
3) r-plyr=1.8.9
4) r-dplyr=1.1.4
5) r-stringr=1.5.1
6) bioconductor-minfi=1.52.0
7) bioconductor-sesame=1.24.0
8) bioconductor-sesamedata=1.24.0
9) conumee2_2.1.2
10) GEOquery_2.74.0
bioconductor-rtracklayer=1.66.0
GenomeInfoDb_1.42.0 
ComplexHeatmap_2.22.0 
r-circlize=0.4.16
plyranges_1.22.0
biomaRt_2.62.0
IlluminaHumanMethylation450kmanifest_0.4.0
IlluminaHumanMethylationEPICmanifest_0.3.0
IlluminaHumanMethylationEPICv2manifest_1.0.0
IlluminaHumanMethylationEPICv2anno.20a1.hg38_1.0.0
IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1 
IlluminaMouseMethylationanno.12.v1.mm10_0.0.2
IlluminaMouseMethylationmanifest_0.0.1

## Input SampleSheet
