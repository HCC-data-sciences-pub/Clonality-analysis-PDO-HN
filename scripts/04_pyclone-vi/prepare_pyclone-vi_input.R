
rm(list = ls())
gc()

## libs 
if(TRUE) {
  
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  
  library(vcfR)
  
}

## ------------------------------------------------

## settings 
if(TRUE) {
  
  project = 'PDO_sm61'
  
  patients.tumor.pdo = c("HN23-10730", "HN24-10810", 
                         "HN24-10856", "HN24-10900", "HN24-10909", "HN24-10936")
  
  print(patients.tumor.pdo)
  
  ## ------------------------------------------------
  
  # my.altDP.cutoff = 0
  # my.altDP.cutoff = 3
  my.altDP.cutoff = 8
}

## ------------------------------------------------

path = 'D:/Users/baor/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61'
setwd(path)

## ------------------------------------------------

## data files 
if(TRUE) {
  
  clinical.file = '../../../sampleinfo/PDO_sm61.exome.sampleinfo.20250102.xlsx'

  ## all samples merged VCF intersect with facets copy number segmentation
  knownVCFallCN.file = 'PDO_sm61.exome.facetsSuite.srt.bed.wknownVCF.bed.slim'
  
  ## tumor purity and copy number seg files
  facetsSuite.file = 'PDO_sm61.exome.facetsSuite.txt'

  ## PASS SNPS: allele specific counts
  # gatkASE.file
  
  ## filtered MAF
  # sommut.file = 'PDO_sm61.exome.fixed.gnomAD_exome_AF0.0001_altDP8.maf'
  sommut.file = 'PDO_sm61.exome.fixed.gnomAD_exome_AF0.0001.maf'
  
}

## ------------------------------------------------

## import and process data 
if(TRUE) {
  
  ## sample info and tumor normal list
  if(TRUE) {
    sampleinfo = data.frame(readxl::read_excel(clinical.file, 
                                               sheet = 'sampleinfo.reorder'
    ),
    stringsAsFactors = F)
    tumornormallist =  data.frame(readxl::read_excel(clinical.file, 
                                                     sheet = 'tumor_normal.reorder'
    ),
    stringsAsFactors = F)
    
    ## ------------------------------------------------
    
    dim(sampleinfo) ## 61 samples total
    dim(tumornormallist) ## 52 tumors 
  }

  ## ------------------------------------------------
  
  ## tumor purity and segs and copy numbers
  if(TRUE) {
    
    facetsSuite = read.delim(facetsSuite.file)
    
    facetsSuite$name = paste0(
      facetsSuite$Tumor,':',facetsSuite$chrom,":",
      facetsSuite$start,':',facetsSuite$end
    )
    
  }
  
  ## ------------------------------------------------
  
  ## copy number of known sites: All
  if(TRUE) {
    knownVCFallCN = read.delim(knownVCFallCN.file,
                               header = F)
    head(knownVCFallCN)
    
    colnames(knownVCFallCN) = c('CHROM','POS','ID','REF','ALT',
                                'QUAL','FILTER',
                                'chrom','chromStart','chromEnd',
                                'name','score','strand','overlap')
    
    knownVCFallCN$Tumor_Sample_Barcode = gsub('[:]\\S+$','',knownVCFallCN$name)
    table(knownVCFallCN$Tumor_Sample_Barcode)
    table(knownVCFallCN$overlap)
    # 1 
    # 636515
    
    ## ------------------------------------------------
    
    dim(knownVCFallCN) ## 636954     15
    knownVCFallCN = knownVCFallCN[knownVCFallCN$Tumor_Sample_Barcode %in%
                                    tumornormallist$Tumor,,drop=F]
    dim(knownVCFallCN) ## 636515     15
    
    knownVCFallCN$Key = paste0(
      knownVCFallCN$CHROM,':',knownVCFallCN$POS,":",
      knownVCFallCN$REF,':',knownVCFallCN$ALT
    )
    
    write.csv(knownVCFallCN,
              file = paste0(project,'.knownVCFallCN',nrow(knownVCFallCN),'.csv'))
  }

  ## ------------------------------------------------
  
  ## patients, gatk ASE counts per sample, merged VCF on PASS SNPs
  if(TRUE) {
    patient.list = list()
    
    ## PASS SNPs: allele specific counts
    gatkASE.list = list()
    
    ## PASS SNPs: per patient merged VCF
    knownVCF.list = list()
    
    ## ------------------------------------------------
    
    for(my.pt in sort(unique(tumornormallist$PID))) {
      
      print(my.pt)
      
      my.df = NULL 
      my.df = tumornormallist[tumornormallist$PID == my.pt,,drop=F]
      print(dim(my.df))
      
      ## keep in mind that I already manually sorted those samples by each patient and by passage number ahead!! no need to resort here. just take the row order as is.
      patient.list[[my.pt]] = as.character(my.df$Tumor)
      
      ## ------------------------------------------------
      
      ## read in merged VCF per pt:
      my.file = NULL 
      my.file = file.path('mutect2_vcf_merged',
                          paste0(my.pt,
                                 '.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz'))
      my.df.2 = NULL
      my.df.2 = read.vcfR(my.file)
      
      class(my.df.2)
      
      length(my.df.2@meta) ## 2859
      dim(my.df.2@fix) ## 897   8
      dim(my.df.2@gt) ## 897  11
      
      knownVCF.list[[my.pt]] = data.frame(my.df.2@fix,
                                          stringsAsFactors = F)
      
      knownVCF.list[[my.pt]]$Key = paste0(
        knownVCF.list[[my.pt]]$CHROM,':',knownVCF.list[[my.pt]]$POS,":",
        knownVCF.list[[my.pt]]$REF,':',knownVCF.list[[my.pt]]$ALT
      )
      
      ## ------------------------------------------------
      
      ## read in gatk ASE counts per sample:
      for(my.sm in my.df$Tumor) {
        print(my.sm)
        
        for(my.tag in c('Pt','All')) {
          
          my.file = NULL 
          my.file = file.path('gatkASEReadCounter',my.tag,
                              paste0(my.sm,
                                     '.bwamem.merged.rmdup.bam.knownSites',
                                     my.tag,'.txt'))
          
          my.df.2 = NULL 
          my.df.2 = read.delim(my.file)
          dim(my.df.2) ## 897  13
          
          my.df.2$Key = paste0(
            my.df.2$contig,':',my.df.2$position,":",
            my.df.2$refAllele,':',my.df.2$altAllele
          )
          
          gatkASE.list[[my.sm]][[my.tag]] = my.df.2
          
        }

      }
      
    }
    
    ## ------------------------------------------------
    
    ## also add the known vcf merged across all tumors 
    knownVCF.list[['All']] = data.frame(read.vcfR(
      'mutect2_vcf_merged/PDO_sm61.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz'
    )@fix,stringsAsFactors = F)
    
    knownVCF.list[['All']]$Key = paste0(
      knownVCF.list[['All']]$CHROM,':',knownVCF.list[['All']]$POS,":",
      knownVCF.list[['All']]$REF,':',knownVCF.list[['All']]$ALT
    )
    
    ## ------------------------------------------------
    
    print(patient.list)
    # $`HN23-10730`
    # [1] "TIIL_0165" "TIIL_0167" "TIIL_0168" "TIIL_0190" "TIIL_0260"
    # 
    # $`HN23-10752`
    # [1] "TIIL_0169" "TIIL_0171" "TIIL_0172" "TIIL_0173" "TIIL_0174" "TIIL_0176"
    # [7] "TIIL_0191" "TIIL_0192"
    # 
    # $`HN24-10791`
    # [1] "TIIL_0177" "TIIL_0179"
    # 
    # $`HN24-10810`
    # [1] "TIIL_0180" "TIIL_0182" "TIIL_0183" "TIIL_0184" "TIIL_0185" "TIIL_0193"
    # [7] "TIIL_0261"
    # 
    # $`HN24-10831`
    # [1] "TIIL_0186" "TIIL_0188" "TIIL_0189" "TIIL_0194" "TIIL_0195"
    # 
    # $`HN24-10856`
    # [1] "TIIL_0262" "TIIL_0264" "TIIL_0265" "TIIL_0266" "TIIL_0267" "TIIL_0268"
    # [7] "TIIL_0269" "TIIL_0270"
    # 
    # $`HN24-10900`
    # [1] "TIIL_0280" "TIIL_0282" "TIIL_0283" "TIIL_0284"
    # 
    # $`HN24-10909`
    # [1] "TIIL_0271" "TIIL_0273" "TIIL_0274" "TIIL_0275" "TIIL_0276" "TIIL_0277"
    # [7] "TIIL_0278" "TIIL_0279"
    # 
    # $`HN24-10936`
    # [1] "TIIL_0285" "TIIL_0287" "TIIL_0288" "TIIL_0289" "TIIL_0290"
    
    print(names(gatkASE.list))
    # [1] "TIIL_0165" "TIIL_0167" "TIIL_0168" "TIIL_0190" "TIIL_0260" "TIIL_0169"
    # [7] "TIIL_0171" "TIIL_0172" "TIIL_0173" "TIIL_0174" "TIIL_0176" "TIIL_0191"
    # [13] "TIIL_0192" "TIIL_0177" "TIIL_0179" "TIIL_0180" "TIIL_0182" "TIIL_0183"
    # [19] "TIIL_0184" "TIIL_0185" "TIIL_0193" "TIIL_0261" "TIIL_0186" "TIIL_0188"
    # [25] "TIIL_0189" "TIIL_0194" "TIIL_0195" "TIIL_0262" "TIIL_0264" "TIIL_0265"
    # [31] "TIIL_0266" "TIIL_0267" "TIIL_0268" "TIIL_0269" "TIIL_0270" "TIIL_0280"
    # [37] "TIIL_0282" "TIIL_0283" "TIIL_0284" "TIIL_0271" "TIIL_0273" "TIIL_0274"
    # [43] "TIIL_0275" "TIIL_0276" "TIIL_0277" "TIIL_0278" "TIIL_0279" "TIIL_0285"
    # [49] "TIIL_0287" "TIIL_0288" "TIIL_0289" "TIIL_0290"
    
    print(names(knownVCF.list))
    # [1] "HN23-10730" "HN23-10752" "HN24-10791" "HN24-10810" "HN24-10831"
    # [6] "HN24-10856" "HN24-10900" "HN24-10909" "HN24-10936" "All"
    
    ## ------------------------------------------------
    
    saveRDS(patient.list,
            file = paste0(project,'.patient.list.rds'))
    saveRDS(gatkASE.list,
            file = paste0(project,'.gatkASE.list.rds'))
    saveRDS(knownVCF.list,
            file = paste0(project,'.knownVCF.list.rds'))
    
    
  }
 
  ## ------------------------------------------------
  
  ## convert sommut to a list of patient => filtered sommmut (from mafm altDP8)
  if(TRUE) {
    sommut = read.delim(sommut.file)
    dim(sommut) ## 11700   222
    data.frame(table(sommut$Tumor_Sample_Barcode)) ## 52 tumors 
    # Var1 Freq
    # 1  TIIL_0165  213
    # 2  TIIL_0167  165
    # 3  TIIL_0168  190
    # 4  TIIL_0169   93
    # 5  TIIL_0171    4
    # .... more!
    
    data.frame(table(sommut$Variant_Classification))
    # Var1 Freq
    # 1                    3'UTR  338
    # 2                  5'Flank   76
    # 3                    5'UTR  311
    # 4      COULD_NOT_DETERMINE    8
    # 5  DE_NOVO_START_OUT_FRAME    7
    # 6          Frame_Shift_Del  181
    # 7          Frame_Shift_Ins   79
    # 8                      IGR   23
    # 9             In_Frame_Del   49
    # 10            In_Frame_Ins    2
    # 11                  Intron 3936
    # 12       Missense_Mutation 4136
    # 13       Nonsense_Mutation  337
    # 14        Nonstop_Mutation    7
    # 15                     RNA  189
    # 16                  Silent 1699
    # 17             Splice_Site  322
    
    data.frame(table(sommut$Variant_Type))
    # Var1  Freq
    # 1  DEL   619
    # 2  DNP   107
    # 3  INS   138
    # 4  SNP 10836
    
    ## ------------------------------------------------
    
    ## sommmut after I manually filtered downstream from MAF files (altDP8)
    print(sommut.file)
    sommut.list = list()
    
    for(my.pt in names(patient.list)) {
      print(my.pt)
      
      my.tumors = NULL 
      my.tumors = patient.list[[my.pt]]
      print(my.tumors)
      
      my.df = NULL 
      my.df = sommut[sommut$Tumor_Sample_Barcode %in% my.tumors,,drop=F]
      dim(my.df) ## 947 222
      
      data.frame(table(my.df$Variant_Type))
      # Var1 Freq
      # 1  DEL   88
      # 2  DNP    5
      # 3  INS   14
      # 4  SNP  840
      
      my.df$Key = paste0(
        my.df$Chromosome,':',my.df$Start_Position,":",
        my.df$Reference_Allele,':',my.df$Tumor_Seq_Allele2
      )
      
      sommut.list[[my.pt]] = my.df
      
    }
    
    print(names(sommut.list))
    # [1] "HN23-10730" "HN23-10752" "HN24-10791" "HN24-10810" "HN24-10831" "HN24-10856"
    # [7] "HN24-10900" "HN24-10909" "HN24-10936"
    
    saveRDS(sommut.list,
            file = paste0(project,'.sommut.list.rds'))
  }
  

}

## ------------------------------------------------

## format variants to pyclone-vi format  
if(TRUE) {
  
  gatkASE.format.list = list()
  pyclone.list = list()
  
  ## ------------------------------------------------
  
  for(my.pt in names(patient.list)) {
    
    print(my.pt)
    ##  "HN23-10730"
    
    my.tumors = NULL 
    my.tumors = patient.list[[my.pt]]
    print(my.tumors)
    
    ## ------------------------------------------------
    
    for(my.tag in c('Pt','All')) {
        
      pyclone.list[[my.pt]][[my.tag]] = NULL 
      
      for(my.sm in my.tumors) {
        
        print(my.sm)
       ##  "TIIL_0165"
        
        ## use pass SNPs
        my.df = NULL 
        my.df = gatkASE.list[[my.sm]][[my.tag]]
        dim(my.df) ## 897  14
        head(my.df)
        
        data.frame(table(my.df$contig)) 
        
        ## exclude sex chrom! hard to interpret downstream
        my.df = my.df[!my.df$contig %in% c('chrX','chrY','chrM'),,drop=F]
        dim(my.df) ## 869  14
        
        # table(my.df$altCount)
        table(my.df$altCount[! my.df$Key %in% sommut.list[[my.pt]]$Key])
        ## lots variants excluded due to low altDP (651 excluded)
        # 0   1   2   3   4   6  12  15  19  24 
        # 506  64  56  16   1   3   1   2   1   1
        ## NULL after I switched to import the maf file with all variants, as expected
        
        ## only keep those kept in downstream analysis (variants passing my filters in any samples from the same patient after my manual filtering, e.g. altDP>=8)
        ## 1/22 update: refactor to be more flexible on altDP cutoffs...
        if(TRUE) {
          
          # my.altDP.cutoff = 0
          # my.altDP.cutoff = 3
          # my.altDP.cutoff = 8
          
          print(my.altDP.cutoff) ## 3

          
          my.select = NULL 
          my.select = sommut.list[[my.pt]]
          dim(my.select) ## 1856  223
          
          my.select = my.select[
              (my.select$t_alt_count +
                 my.select$t_ref_count) >= 10 &
              (my.select$n_alt_count +
                 my.select$n_ref_count) >= 10,,drop=F]
          dim(my.select) ## 1687  223
          
          my.select = my.select[my.select$t_alt_count >= my.altDP.cutoff,,drop=F]
          dim(my.select) ## 1610  223

          my.df = my.df[my.df$Key %in% my.select$Key,,drop=F]
          dim(my.df) ## altDP cutpff=8: 216 ; cutoff=0: 766; cutoff=3: 705
          
          
        }
        
        ## ------------------------------------------------
        
        ## gene name and annotation!!
        
        my.anno = sommut.list[[my.pt]]
        head(my.anno$Key)
        
        my.anno = unique(my.anno[,c("Key",
                                    'Hugo_Symbol',
                                    'Variant_Classification',
                                    'Variant_Type',
                                    'Annotation_Transcript',
                                    'cDNA_Change',
                                    'Protein_Change')])
        my.anno[my.anno==''] = '.'
        
        dim(my.anno) ## 1141    7
        sum(duplicated(my.anno$Key)) ## 0; all uniq!!
        
        dim(my.df) ## 705  14
        my.df = merge(my.df, my.anno, by = 'Key')
        dim(my.df) ## 705  20
        head(my.df)

        ## ------------------------------------------------
        
        ## add copy number 
        
        my.facets = NULL 
        my.facets = facetsSuite[facetsSuite$Tumor==my.sm,,drop=F]
        dim(my.facets) ## 60 24
        head(my.facets)
        
        my.df.2 = NULL 
        my.df.2 = knownVCFallCN[knownVCFallCN$Tumor_Sample_Barcode == my.sm,,drop=F]
        dim(my.df.2) ## 12236    16
        head(my.df.2)
        
        sum(duplicated(my.facets$name)) ## 0; must be unique!!
        my.facets = merge(my.facets,my.df.2[,c('Key','name','Tumor_Sample_Barcode')],
                          by = 'name')
        dim(my.facets) ## 12236    26
        head(my.facets)
        sum(duplicated(my.facets$name)) ## 12179; dup exists ,as expected
        
        ## test purpose only... why 2 mutations missing afer merging???
        if(FALSE) {
          
          w = NULL 
          w = merge(my.df, my.facets, by = 'Key')
          
          my.df[!my.df$Key %in% w$Key,]
          # contig  position variantID refAllele altAllele refCount altCount totalCount
          # 529  chr11 134387630         .         G         A      244        1        245
          # 659  chr16     16893         .         C         T      143        0        143
          # lowMAPQDepth lowBaseQDepth rawDepth otherBases improperPairs
          # 529            0             1      253          0             7
          # 659            0             0      145          0             2
          # Key
          # 529 chr11:134387630:G:A
          # 659     chr16:16893:C:T
          
          ## still do not know why those variants were not present in the copy number file :(
        }
        
        dim(my.df) ## 869  14; 216  14 ; 705  20
        head(my.df)
        sum(duplicated(my.df$Key)) ## 0; must be unique!!
        my.df = merge(my.df, my.facets, by = 'Key')
        dim(my.df) ## 867 39; 215  39; 703  45
        
        head(my.df)
        sum(duplicated(my.df$Key)) ## 0; still unique, as expected
        
        ## ------------------------------------------------
        
        ## add tumor purity (or content)
        
        my.purity = NULL 
        my.purity = unique(facetsSuite$purity[facetsSuite$Tumor==my.sm])
        my.purity = as.numeric(as.character(my.purity))
        print(my.purity) ## 0.309253
        
        my.df$purity = my.purity
        
        ## ------------------------------------------------
        
        ## exclude variants when major < minor cn or cn is NA
        dim(my.df)
        my.df = my.df[!is.na(my.df$tcn.em) &
                        !is.na(my.df$lcn.em),,drop=F]
        dim(my.df) ## 703  45
        
        x = NULL 
        x = which((my.df$tcn.em - my.df$lcn.em) >= my.df$lcn.em)
        length(x)
        
        my.df = my.df[x,,drop=F]
        dim(my.df) ## 703  45
        
        ## ------------------------------------------------
        
        ## format to pyclone input cols 
        
        my.df.format = NULL 
        my.df.format = 
          data.frame(
            
            mutation_id = paste0(my.df$Key,"!",
                                 my.df$Hugo_Symbol,":",
                                 my.df$Protein_Change,":",
                                 my.df$Variant_Classification,":",
                                 my.df$Variant_Type),
            sample_id = my.sm,
            ref_counts = my.df$refCount,
            alt_counts = my.df$altCount,
            major_cn = my.df$tcn.em - my.df$lcn.em,
            minor_cn = my.df$lcn.em,
            normal_cn = 2,
            tumour_content = my.df$purity,
            error_rate = 0.001, ## default value
            
            stringsAsFactors = F
            
          )
        
        dim(my.df.format) ## 867   9
        head(my.df.format)
        
        ## ------------------------------------------------
        
        gatkASE.format.list[[my.sm]][[my.tag]] = my.df
        
        pyclone.list[[my.pt]][[my.tag]] = rbind(
          pyclone.list[[my.pt]][[my.tag]],
          my.df.format
        )
      }
      
      
      
    }
    
  }
  
  ## ------------------------------------------------
  
  print(names(gatkASE.format.list))
  # [1] "TIIL_0165" "TIIL_0167" "TIIL_0168" "TIIL_0190" "TIIL_0260" "TIIL_0169"
  # [7] "TIIL_0171" "TIIL_0172" "TIIL_0173" "TIIL_0174" "TIIL_0176" "TIIL_0191"
  # [13] "TIIL_0192" "TIIL_0177" "TIIL_0179" "TIIL_0180" "TIIL_0182" "TIIL_0183"
  # [19] "TIIL_0184" "TIIL_0185" "TIIL_0193" "TIIL_0261" "TIIL_0186" "TIIL_0188"
  # [25] "TIIL_0189" "TIIL_0194" "TIIL_0195" "TIIL_0262" "TIIL_0264" "TIIL_0265"
  # [31] "TIIL_0266" "TIIL_0267" "TIIL_0268" "TIIL_0269" "TIIL_0270" "TIIL_0280"
  # [37] "TIIL_0282" "TIIL_0283" "TIIL_0284" "TIIL_0271" "TIIL_0273" "TIIL_0274"
  # [43] "TIIL_0275" "TIIL_0276" "TIIL_0277" "TIIL_0278" "TIIL_0279" "TIIL_0285"
  # [49] "TIIL_0287" "TIIL_0288" "TIIL_0289" "TIIL_0290"
  
  print(names(pyclone.list))
  # [1] "HN23-10730" "HN23-10752" "HN24-10791" "HN24-10810" "HN24-10831" "HN24-10856"
  # [7] "HN24-10900" "HN24-10909" "HN24-10936"
  
  saveRDS(gatkASE.format.list,
          file = paste0(project,'.gatkASE.format.list.rds'))
  saveRDS(pyclone.list,
          file = paste0(project,'.gatkASE.pyclone.list.rds'))
  
}

## ------------------------------------------------

## print pyclone inputs to local files  
if(TRUE) {
  
  my.outdir.0 = NULL 
  my.outdir.0 = paste0('pyclone-vi_altDP',my.altDP.cutoff)
  
  print(my.outdir.0)
  
  if (! dir.exists(my.outdir.0)) { dir.create(my.outdir.0, recursive = T)}
  
  ## ------------------------------------------------
  
  for(my.pt in names(patient.list)) {
    
    print(my.pt)
    
    my.tumors = NULL 
    my.tumors = patient.list[[my.pt]]
    print(my.tumors)
    
    ## ------------------------------------------------
    
    for(my.tag in c('Pt','All')) {
      
      my.outdir = file.path(my.outdir.0,my.tag)
      
      if (! dir.exists(my.outdir)) { dir.create(my.outdir, recursive = T)}
      
      
      my.df = NULL 
      my.df = pyclone.list[[my.pt]][[my.tag]]
      
      write.table(my.df,
                  file = file.path(my.outdir,
                                   paste0(project,'.gatkASE.pyclone.list.',
                                my.pt,'.',my.tag,'.txt')),
                  sep = '\t', quote = F, row.names = F, col.names = T)
      
    }
    
  }
  
}


## ------------------------------------------------

sessionInfo()

