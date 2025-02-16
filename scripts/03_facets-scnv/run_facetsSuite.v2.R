
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
  
  facetsSuite.out = NULL 
  facetsSuite.out.list = list()
  
  # my.sm = 'PNATM-101-tumor'
  for(my.sm in tumors) {
    
    print(my.sm)
    
    ## --------------------------------------------
    
    datafile = NULL 
    rcmat = NULL 
    
    datafile = paste0('facets/',
                      my.sm,'.bwamem.merged.rmdup.bam.snp_pileup.csv.gz')
    rcmat = readSnpMatrix(datafile) ## note that using read_snp_matrix (fuc from facetssuite) does not work in this case. must use readSnpMatrix  (func from facets). because I used snp-mpileup function from origial facets (not facetssuite) to generate the snp gz file
    head(rcmat)
    
    ## --------------------------------------------
    
    my.out = NULL 
    my.out = run_facets(rcmat, cval = 150, genome = 'hg38',seed = 1234)
    
    names(my.out)
    # [1] "snps"      "segs"      "purity"    "ploidy"    "dipLogR"   "alBalLogR"
    # [7] "flags"     "em_flags"  "loglik" 
    
    ## --------------------------------------------
    
    facetsSuite.out = rbind(facetsSuite.out,
                       data.frame(Tumor = my.sm,
                                  purity = my.out$purity,
                                  ploidy = my.out$ploidy,
                                  dipLogR = my.out$dipLogR,
                                  flags = paste(my.out$flags,
                                                collapse = ';'),
                                  em_flags = paste(my.out$em_flags,
                                                collapse = ';'),
                                  my.out$segs,
                                  stringsAsFactors = F))
    
    facetsSuite.out.list[[my.sm]] = my.out
    
    
  }
  
  saveRDS(facetsSuite.out.list,
          file.path(paste0(project,'.',data.type,'.facetsSuite.rds')))
  
  write.csv(facetsSuite.out,
            file.path(paste0(project,'.',data.type,'.facetsSuite.csv')))
  
  write.table(facetsSuite.out,
            file.path(paste0(project,'.',data.type,'.facetsSuite.txt')),
            sep = '\t', quote = F, row.names = F, col.names = T)
  
  ## --------------------------------------------
  
  colnames(facetsSuite.out)
  x = NULL 
  x = unique(facetsSuite.out[,c('Tumor','purity','ploidy','flags')])
  
  write.csv(x,
            file.path(paste0(project,'.',data.type,'.facetsSuite.purity.csv')))
  

}

## --------------------------------------------

## annotating .... only ran once !
if(FALSE) {
  
  facetsSuite.out.list = readRDS(file.path(
    paste0(project,'.',data.type,'.facetsSuite.rds')))
  
  names(facetsSuite.out.list)
  
  maf = read.delim('PDO_sm61.exome.fixed.gnomAD_exome_AF0.0001_altDP8.maf')
  table(maf$Tumor_Sample_Barcode)
  
  ## --------------------------------------------
  
  ## save various facetssuite functions ....
  facetsSuite.func.list = list()

  for(my.sm in sort(names(facetsSuite.out.list))) {
    
    print(my.sm)
    
    my.out = NULL 
    my.out = facetsSuite.out.list[[my.sm]]
    
    ## --------------------------------------------
    
    ## run various facetssuite functions ....
    
    for(my.method in c("em", "cncf")) {
      print(my.method)
      
      ## Arm-level changes:
      ## Get the altered chromosome arms in sample. Does not include the acrocentric p arms of chromosomes 12, 14, 15, 31, and 22. 
      my.out.2 = NULL 
      my.out.2 = arm_level_changes(my.out$segs, my.out$ploidy, genome = "hg38",
                        algorithm = my.method)
      
      facetsSuite.func.list[[my.method]][[my.sm]][['arm_level_changes']] = my.out.2
      
      ## --------------------------------------------
      
      ## Estimate CCFs of somatic mutations
      my.out.2 = NULL 
      my.maf = NULL 
      my.maf = maf[maf$Tumor_Sample_Barcode==my.sm,,drop=F]
      if(nrow(my.maf) > 1) {
        
        my.out.2 = ccf_annotate_maf(my.maf, my.out$segs, my.out$purity, 
                                    algorithm = my.method)
        table(my.out.2$clonality)
      }
      
      facetsSuite.func.list[[my.method]][[my.sm]][['ccf_annotate_maf']] = my.out.2
      
      ## --------------------------------------------
      
      ## Sample QC
      my.out.2 = NULL 
      my.out.2 = check_fit(my.out, genome = "hg38",
                algorithm = my.method)	
      
      my.maf = NULL 
      my.maf = maf[maf$Tumor_Sample_Barcode==my.sm,,drop=F]
      if(nrow(my.maf) > 1) {
        
        my.out.2 = check_fit(my.out, genome = "hg38",
                             algorithm = my.method, maf = my.maf)	
      }
      
      facetsSuite.func.list[[my.method]][[my.sm]][['check_fit']] = my.out.2
      
      ## --------------------------------------------
      
      ## 	Copy-number based scores
      # copy_number_scores
      # Calculate the following
      # 
      # Fraction of genome altered:
      #   
      #   Fraction of genome altered and genome doubling flag.
      # Fraction LOH:
      #   
      #   Fraction of genome with LOH and flag for hypoploidy.
      # LST score:
      #   
      #   Large-scale state transitions, see source URL.
      # NtAI:
      #   
      #   Telomeric allelic imbalance, see source URL.
      # HRD-LOH:
      #   
      #   HRD-LOH score, see source URL.
      
      my.out.2 = NULL 
      my.out.2[['calculate_fraction_cna']] = data.frame(
        calculate_fraction_cna(my.out$segs, my.out$ploidy, genome = "hg38",
                               algorithm = my.method)
      )
      my.out.2[['calculate_loh']] = data.frame(
        calculate_loh(my.out$segs, my.out$snps, genome = "hg38",
                      algorithm = my.method, hypoploidy_threshold = 0.5)
      )
      my.out.2[['calculate_ntai']] = data.frame(
        calculate_ntai(my.out$segs, my.out$ploidy, genome = "hg38",
                       algorithm = my.method, min_size = 0, min_probes = 250)
      )
      my.out.2[['calculate_lst']] = data.frame(
        calculate_lst(my.out$segs, my.out$ploidy, genome = "hg38",
                      algorithm = my.method, min_size = 1e+07, min_probes = 50)
      )
      my.out.2[['calculate_fraction_cna']] = data.frame(
        calculate_hrdloh(my.out$segs, my.out$ploidy, algorithm = my.method)
      )
      
      my.out.3 = NULL 
      for(my.func in names(my.out.2)) {
        colnames(my.out.2[[my.func]]) = 
          paste0(my.func,'.',colnames(my.out.2[[my.func]]))
        
        my.out.3 = rbind(my.out.3,
                         data.frame(t(my.out.2[[my.func]])))
      }
      colnames(my.out.3)[1] = 'value'
      
      facetsSuite.func.list[[my.method]][[my.sm]][['copy_number_scores']] = my.out.3

      ## --------------------------------------------
      
      ## Format '.seg' file
      # format_igv_seg	
      my.out.2 = NULL 
      my.out.2 = format_igv_seg(my.out, my.sm, normalize = TRUE)
      facetsSuite.func.list[[my.method]][[my.sm]][['format_igv_seg']] = my.out.2
      
      ## --------------------------------------------
      
      ## Get gene-level changes in copy number from FACETS output.
      # gene_level_changes	
      my.out.2 = NULL 
      my.out.2 = gene_level_changes(my.out, genome =  "hg38",
                                    algorithm = my.method)
      dim( my.out.2) ## 8166   24
      facetsSuite.func.list[[my.method]][[my.sm]][['gene_level_changes']] = my.out.2

      ## --------------------------------------------

      ## Plot FACETS output
      # plot_facets	
#         Log-ratio: Copy-number segmentation based on tumor-normal read coverage comparison.
#       
#       Allelic imbalance: Allelic imbalance based on somatic changes in zygosity of heterozygous SNPs in the normal.
#       
#       Integer copy number: Inference of the major and minor copy-number states based on the tumor-normal log ratio and allelic content.
#       
#       Cellular fraction: Estimate of fraction of cells in sample harboring each segment.
  
      p.multi = list()
      p.multi[['cnlr_plot']] = cnlr_plot(my.out, genome = "hg38")
      p.multi[['valor_plot']] = valor_plot(my.out, genome = "hg38")
      p.multi[['cf_plot']] = cf_plot(my.out, genome = "hg38")
      p.multi[['icn_plot']] = icn_plot(my.out, genome = "hg38")
      
      # p.multi[['closeup_plot']] = closeup_plot(my.out, genome = "hg38")
      ## Error: Specify at least a gene or chromosome to zoom in on.
      p.multi[['closeup_plot']] = closeup_plot(my.out, genome = "hg38", plot_chroms = c('6'))

      pdf(file = paste0(project,'.',data.type,'.facetsSuite.',
                        my.method,'.',my.sm,'.plot.pdf'),
          height = 10,width = 10)
      print(p.multi[[1]] + p.multi[[2]] + p.multi[[3]] +
              p.multi[[4]] + 
              plot_layout(ncol = 1))
      dev.off()
      
      pdf(file = paste0(project,'.',data.type,'.facetsSuite.',
                        my.method,'.',my.sm,'.closeup_plot.chr6.pdf'),
          height = 5,width = 10)
      print(p.multi[['closeup_plot']])
      dev.off()

      
      ## --------------------------------------------
      
      # ## Plot segmentation profile
      # # plot_segmentation	
      # # Creates an IGV-like graphical representation of the copy-number segments across the samples in a segmentation object, as output by run_facets.
      # my.segs.igv = NULL 
      # my.segs.igv = format_igv_seg(my.out, my.sm, normalize = TRUE)
      # plot_segmentation(my.segs.igv, 
      #                   sample_order = c(my.sm),
      #                   cap_log_ratios = TRUE)
      # # Error in `mutate()`:
      # #   i In argument: `.data`.
      # # Caused by error:
      # #   ! `.data` must be a vector, not a <rlang_data_pronoun> object.
      
      ## --------------------------------------------
      
      ## other funcs ; not useful here
      # ## Read counts file
      # read_snp_matrix	
      # ## 	Run FACETS
      # run_facets
      # ## 	Test FACETS output
      # test_facets_output
      # ## 	Test MAF
      # test_maf
      # ## 	Test read counts
      # test_read_counts
      
    }
    
    
    


    
    
    
  }
  
  saveRDS(facetsSuite.func.list,
          file.path(paste0(project,'.',data.type,'.facetsSuite.func.rds')))
  
}



## --------------------------------------------













