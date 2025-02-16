
rm(list = ls())
gc()

## -----------------------------------------------

## libs 
if(TRUE) {
  # devtools::install_github("mskcc/facets", build_vignettes = F)
  # devtools::install_github("mskcc/facets-suite")
  
  library(facets)
  library(facetsSuite)
}

## -----------------------------------------------

## settings
if(TRUE) {
  
  project = 'PDO_sm61'
  
  data.type = 'exome'
  
}

## -----------------------------------------------

path = file.path('D:/Users/baor/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/',data.type,project)
setwd(path)

## -----------------------------------------------

## data files 
if(TRUE) {
  
  tumornormallist.file = ''
  
  if(data.type =='exome') {
    tumornormallist.file = '../../../sampleinfo/PDO_sm61.exome.tumor_normal.20250102.list'
  } 
  # if(data.type =='genome') {
  #   tumornormallist.file = '../sampleinfo/ATM_PAAD.genome.tumor_normal.20200526.list'
  # }
  
}

## -----------------------------------------------

## import data 
if(TRUE) {
  
  tumornormallist = read.delim(tumornormallist.file, stringsAsFactors = F)
  tumors = sort(unique(tumornormallist$Tumor))
  
  print(tumors)
  
}

## --------------------------------------------

## only ran once !
if(FALSE) {
  
  facets.out = NULL 
  facets.out.list = list()
  
  # my.sm = 'PNATM-101-tumor'
  for(my.sm in tumors) {
    
    print(my.sm)
    
    set.seed(1234)
    rm(datafile, rcmat, xx, oo, fit, p1, p2)
    
    
    datafile = paste0('facets/',
                      my.sm,'.bwamem.merged.rmdup.bam.snp_pileup.csv.gz')
    rcmat = readSnpMatrix(datafile)
    xx = preProcSample(rcmat,gbuild = 'hg38')
    
    oo=procSample(xx,cval=150)
    
    oo$dipLogR
    
    fit=emcncf(oo)
    
    head(fit$cncf)
    
    facets.out = rbind(facets.out,
                       data.frame(Tumor = my.sm,
                                  purity = fit$purity,
                                  ploidy = fit$ploidy,
                                  dipLogR = fit$dipLogR,
                                  # emflags = fit$emflags,
                                  fit$cncf,
                                  stringsAsFactors = F))
    
    facets.out.list[[my.sm]] = fit
    
    
    fit$purity
    fit$ploidy
    
    pdf(file = file.path(paste0(my.sm,'.facets.fit.p1.pdf')), height = 6, width = 10)
    plotSample(x=oo,emfit=fit)
    dev.off()
    
    pdf(file = file.path(paste0(my.sm,'.facets.fit.p2.pdf')), height = 6, width = 10)
    logRlogORspider(oo$out, oo$dipLogR)
    dev.off()
    
  }
  
  saveRDS(facets.out.list,
          file.path(paste0(project,'.',data.type,'.facets.fit.rds')))
  
  write.csv(facets.out,
            file.path(paste0(project,'.',data.type,'.facets.fit.csv')))
  
  write.table(facets.out,
              file.path(paste0(project,'.',data.type,'.facets.fit.txt')),
              sep = '\t', quote = F, row.names = F, col.names = T)
  
  ## --------------------------------------------
  
  colnames(facets.out)
  x = NULL 
  x = unique(facets.out[,c('Tumor','purity','ploidy')])
  
  write.csv(x,
            file.path(paste0(project,'.',data.type,'.facets.fit.purity.csv')))
  
}

## --------------------------------------------

## annotating ....
if(TRUE) {
  
  facets.out = read.csv(file.path(paste0(project,'.',data.type,'.facets.fit.csv')),
                        row.names = 1, stringsAsFactors = F)
  
  table(facets.out$ploidy)
  
  table(facets.out$Tumor)
  
  ## classify LOH events etc.
  ## https://github.com/mskcc/facets/issues/62
  # tcn.em	lcn.em	type
  # 0	0	DEL
  # 0	NA	DEL
  # 1	0	HEMIZYGOTE
  # 1	NA	HEMIZYGOTE
  # 2	0	LOH
  # 2	1	NEUTR
  # 2	NA	NEUTR/Unknown
  # 2+	0	DUP-LOH
  # 2+	1+	DUP
  # 2+	NA	DUP-[LOH?]
  
  sum(is.na(facets.out$tcn.em))
  sum(is.na(facets.out$lcn.em)) ## 627
  
  if(TRUE) {
    
    facets.out$type = 'Unknown'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            !is.na(facets.out$lcn.em) &
                            facets.out$tcn.em == 0 &
                            facets.out$lcn.em == 0)] = 'DEL'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            is.na(facets.out$lcn.em) &
                            facets.out$tcn.em == 0)] = 'DEL'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            !is.na(facets.out$lcn.em) &
                            facets.out$tcn.em == 1 &
                            facets.out$lcn.em == 0)] = 'HEMIZYGOTE'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            is.na(facets.out$lcn.em) &
                            facets.out$tcn.em == 1)] = 'HEMIZYGOTE'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            !is.na(facets.out$lcn.em) &
                            facets.out$tcn.em == 2 &
                            facets.out$lcn.em == 0)] = 'LOH'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            !is.na(facets.out$lcn.em) &
                            facets.out$tcn.em == 2 &
                            facets.out$lcn.em == 1)] = 'NEUTR'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            is.na(facets.out$lcn.em) &
                            facets.out$tcn.em == 2)] = 'NEUTRorUnknown'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            !is.na(facets.out$lcn.em) &
                            facets.out$tcn.em ==2 &
                            facets.out$lcn.em == 0)] = 'LOH'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            !is.na(facets.out$lcn.em) &
                            facets.out$tcn.em >2 &
                            facets.out$lcn.em == 0)] = 'DUP_LOH'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            !is.na(facets.out$lcn.em) &
                            facets.out$tcn.em >2 &
                            facets.out$lcn.em >1)] = 'DUP'
    facets.out$type[which(!is.na(facets.out$tcn.em) &
                            is.na(facets.out$lcn.em) &
                            facets.out$tcn.em >2 )] = 'DUP_LOH_maybe'
    
    table(facets.out$type)
    # DEL            DUP        DUP_LOH  DUP_LOH_maybe     HEMIZYGOTE 
    # 74            578            167            390            340 
    # LOH          NEUTR NEUTRorUnknown        Unknown 
    # 394            798            177            461 
    
    
  }
  
  x = NULL 
  x = grep('LOH',facets.out$type)
  
  write.csv(facets.out,
            file = file.path(paste0(project,'.',data.type,'.facets.fit.wType.csv')))
  
  write.table(facets.out,
            file = file.path(paste0(project,'.',data.type,'.facets.fit.wType.txt')),
            sep = '\t', quote = F, row.names = F, col.names = T)
  
  write.csv(facets.out[x,],
            file = file.path(paste0(project,'.',data.type,'.facets.fit.wType.LOH.csv')))
  
  write.table(facets.out[x,],
              file = file.path(paste0(project,'.',data.type,'.facets.fit.wType.LOH.txt')),
              sep = '\t', quote = F, row.names = F, col.names = T)
  
  sum(is.na(facets.out$type))
  table(facets.out$type)
  # DEL           DUP       DUP_LOH DUPorLOHmaybe    HEMIZYGOTE 
  # 74          1837           561           567           340 
  
  ## --------------------------------------------
  
  if(FALSE) {
    ## ATM 
    # chr11   108222831       108369102       ATM     .       +
    atm.chr = '11'
    atm.start = 108222831
    atm.end = 108369102
    facets.out.atm = facets.out[
      
      facets.out$chrom == atm.chr &
        (
          (
            facets.out$start <= atm.start &
              facets.out$end >= atm.end 
          ) |
            (
              facets.out$start <= atm.start &
                facets.out$end>= atm.start &
                facets.out$end <= atm.end
            ) |
            (
              facets.out$start >= atm.start &
                facets.out$start <= atm.end &
                facets.out$end >= atm.end
            ) |
            (
              facets.out$start >= atm.start &
                facets.out$end <= atm.end 
            )
        )
      ,
    ]
    
    dim(facets.out.atm)
    
    write.csv(facets.out.atm,
              file = file.path(paste0(project,'.',data.type,'.facets.fit.wType.atm.csv')))
  }
  
  
}

# 
# 
# ## --------------------------------------------
# 
# facets.out.list = readRDS(file.path(
#                                     paste0(project,'.',data.type,'.facets.fit.rds')))
# 
# for(my.sm in sort(names(facets.out.list))) {
# 
#   print(my.sm)
#   
#   fit = NULL 
#   fit = facets.out.list[[my.sm]]
# 
#   my.df = fit$cncf
# 
#   ## gene centric annotation
#   my.df.1 = gene_level_changes(fit, genome = 'hg38',
#                                algorithm = 'em')
#   my.df.2 = gene_level_changes(fit, genome = 'hg38',
#                                algorithm = 'cncf')
# 
# 
# 
# }









