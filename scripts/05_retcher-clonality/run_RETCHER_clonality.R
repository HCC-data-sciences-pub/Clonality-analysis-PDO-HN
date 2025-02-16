
## on htc.... on command line 

# cd /ihome/rbao/rib37 ## go to home directory
# module load gcc/12.2.0 r/4.3.0
# R ## open console on command line 

## within R console, run following statements 

## ----------------------------------------

## within console, run following statements 

rm(list = ls())
gc()

## --------------------------------------------

## https://stackoverflow.com/questions/28633600/r-rgl-package-error-while-loading-the-library
options(rgl.useNULL=TRUE)

# The rgl.useNULL = TRUE option in R is used to specify how the rgl package (used for 3D visualizations) handles the rendering of 3D plots. Specifically, setting rgl.useNULL = TRUE ensures that rgl does not attempt to open an interactive 3D rendering window, which is particularly useful in non-interactive or headless environments

## --------------------------------------------

## https://github.com/zlsys3/RETCHER

# devtools::install_github("zlsys3/RETCHER")

library(RETCHER)
# install_sciClone() ## I already installed it before 

## --------------------------------------------

## take in command line arguments ...
if(FALSE) {
  
  args = commandArgs(trailingOnly=TRUE)
  
  if (length(args)==0) {
    stop("At least one argument must be supplied (my.pt).n", call.=FALSE)
  } else if (length(args)>0) {
    
    my.pt = as.character(args[1])
    
  
  }
  
}

## --------------------------------------------

## settings
if(TRUE) {
  
  # outdir = paste0("./RETCHER_out_",my.pt,'_17')
  outdir = paste0("./RETCHER_out_",my.pt)
  
  print(my.pt)
  print(outdir)
  
  ## --------------------------------------------
  
  altDP.cutoff = 3 ## 0 3 8
  
  ## --------------------------------------------
  
  patients.tumor.pdo = c("HN23-10730", "HN24-10810", 
                         "HN24-10856", "HN24-10900", "HN24-10909", "HN24-10936")
  
  print(patients.tumor.pdo)
  
  ## --------------------------------------------
  
  sample.ex = c('TIIL_0282')
  
  print(sample.ex)
  
}

## --------------------------------------------

## input format 
# columns         description
# chr     The chromosome name of the mutation
# pos     The position of the mutation
# gene    The identifier of the mutation
# ref     The number of reads for the reference allele of the mutation
# alt     The number of reads for the alternate allele of the mutation
# major   The major copy number of the mutation
# minor   The minor copy number of the mutation
# purity  The purity of the sample

## test only; example run
if(FALSE) {
  
  data(simdata)
  head(simdata$sim_s1)
  
  pipeline_res <- RunPipeline(inputList = simdata)
  # or
  pipeline_res <- RunPipeline(inputList = simdata, saveDir = "./results")
  
}

## --------------------------------------------

## altDP3, altDP8, altDP0
path = paste0('C:/Users/Owner/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61/pyclone-vi_altDP',altDP.cutoff,'/Pt/')
# path = '/ix/jluke/Projects/sunny/TIIL-017-RoseProj/results/exome/PDO_sm61_multisample_grch38/pyclone-vi_altDP3/Pt'
setwd(path)

print(path)

## --------------------------------------------

# ## redirect stdout into local files 
# ## https://stackoverflow.com/questions/7096989/how-to-save-all-console-output-to-file-in-r

# con = NULL 
# con <- file("RETCHER.20250122.log")
# sink(con, append=TRUE)
# sink(con, append=TRUE, type="message")

## --------------------------------------------


# my.pt = patients.tumor.pdo[1]
# print(my.pt)


# for(my.pt in patients.tumor.pdo) {

print(my.pt)

my.file = NULL 
my.file = paste0('PDO_sm61.gatkASE.pyclone.list.',my.pt,'.Pt.txt')

my.df = NULL 
my.df = read.delim(my.file)

my.list = list()
for(my.sm in sort(unique(my.df$sample_id))) {
  
  print(my.sm)
  
  ## --------------------------------------------
  
  if(my.sm %in% sample.ex) {
    
    print("sample needs to be excluded. skip")
    next
    
  }
  
  ## --------------------------------------------
  
  my.df.2 = NULL 
  my.df.2 = my.df[my.df$sample_id==my.sm,,drop=F]
  head(my.df.2)
  
  my.df.2$mutation_id.0 = my.df.2$mutation_id
  
  ## chr1:100720686:C:G!VCAM1:p.S92C:Missense_Mutation:SNP
  my.df.2 = 
    tidyr::separate(my.df.2, col = mutation_id, sep = ':|!',
                    into = c('CHROM','POS','REF','ALT',
                             'Hugo_Symbol',
                             'Protein_Change',
                             'Variant_Classification',
                             'Variant_Type'
                    ))
  
  # columns         description
  # chr     The chromosome name of the mutation
  # pos     The position of the mutation
  # gene    The identifier of the mutation
  # ref     The number of reads for the reference allele of the mutation
  # alt     The number of reads for the alternate allele of the mutation
  # major   The major copy number of the mutation
  # minor   The minor copy number of the mutation
  # purity  The purity of the sample
  
  my.df.2.format = NULL 
  my.df.2.format = 
    data.frame(
      chr = as.character(my.df.2$CHROM),
      pos = as.numeric(as.character(my.df.2$POS)),
      gene = as.character(my.df.2$mutation_id.0),
      ref = as.numeric(as.character(my.df.2$ref_counts)),
      alt = as.numeric(as.character(my.df.2$alt_counts)),
      major = as.numeric(as.character(my.df.2$major_cn)),
      minor = as.numeric(as.character(my.df.2$minor_cn)),
      purity = as.numeric(as.character(my.df.2$tumour_content)),
      
      stringsAsFactors = F
      
    )
  
  head(my.df.2.format)
  
  print(dim(my.df.2.format))
  my.df.2.format = my.df.2.format[!is.na(my.df.2.format$major) &
                                    my.df.2.format$major >0 &
                                    !is.na(my.df.2.format$minor) &
                                    my.df.2.format$minor >=0
                                  ,,drop=F]
  print(dim(my.df.2.format))
  
  my.list[[my.sm]] = my.df.2.format
  
}

## --------------------------------------------

pipeline_res = NULL 
pipeline_res = RunPipeline(inputList = my.list,
                           saveDir = outdir,
                           calcCCFMethod = "distance", ## default = distance
                           maxCCFValue = 1.2,
                           # sampleNames = NULL,
                           minimumDepth = 10, ## default = 20 (total depth)
                           maximumClusters = 10, ## default = 10
                           useSexChrs = FALSE, ## default = TRUE
                           bootstrapIter = 2000,
                           bootstrapConf = 0.95,
                           usePercentage = TRUE,
                           LimitCCF = TRUE,
                           clusterMinSnvNum = 5, ## default = 5
                           removeOutliers = TRUE, ## default = TRUE
                           mergeUndiffClust = TRUE,
                           mergePvalue = 0.05, ## default = 0.05
                           error_buf = 0.1,
                           enumAllTree = TRUE,
                           maxTreeNum = 20,
                           allowRemoveCluster = TRUE, ## default = TRUE
                           minPresenceCCF = 0.01 ## default = 0.05
                           
)


## HN24-10856 errored out when runing with distance:
## Error in if (distance < min.distance) { : 
# missing value where TRUE/FALSE needed
## changed to threshold and it runs 
## then I decided to just switch from distance to threshold, because other cases threw out warnings too with distance?! not sure why

## --------------------------------------------

saveRDS(my.list,
        file = file.path(outdir,
                         paste0('pipeline_res.in.rds')))

saveRDS(pipeline_res,
        file = file.path(outdir,
                         paste0('pipeline_res.out.rds')))

# }

## --------------------------------------------

sessionInfo()

## --------------------------------------------

# # Restore output to console
# sink() 
# sink(type="message")

# ## And look at the log...
# # cat(readLines("RETCHER.20250122.log"), sep="\n")



