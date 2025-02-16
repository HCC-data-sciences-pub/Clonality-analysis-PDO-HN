
rm(list = ls())
gc()

## reduce altDP filter from 8 to 5

if(TRUE) {
  library(cluster)
  library(Biobase)
  library(ctc)
  library(ape)
  library(ggdendro)
  library(ggplot2)
  library(reshape2)
  library(grid)
  library(gplots)
  library(stringr)
  library(RColorBrewer)
  
  library(ComplexHeatmap)
  library(broom)
  library(survival)
  library(survminer)
  library(circlize)
  
  ## https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
  ## https://github.com/cran/ggsci/blob/master/R/discrete-nejm.R
  library(ggsci)
  library(scales)
  library(ggpubr)
  library(plotrix)
  library(maftools)
  library(edgeR)
  library(xCell)
  library(fBasics)
  library(data.table)
}

## --------------------------------------------------------

## settings
if(TRUE) {
  

  sample.ex = c()
  
  output = 'PDO_sm61.exome' 
  
  var.funcs = c('DE_NOVO_START_IN_FRAME','DE_NOVO_START_OUT_FRAME',
                'Frame_Shift_Del','Frame_Shift_Ins',
                'In_Frame_Del','In_Frame_Ins',
                'Missense_Mutation','Nonsense_Mutation',
                'Splice_Site','Translation_Start_Site')
  
}

## --------------------------------------------------------

path = 'D:/Users/baor/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61'
setwd(path)

## --------------------------------------------------------

## only run once: filter variants 
if(FALSE) {

  ## import variants 
  if(TRUE) {
    
    clinical.file = '../../../sampleinfo/PDO_sm61.exome.sampleinfo.20250102.xlsx'
    
    sampleinfo = data.frame(readxl::read_excel(clinical.file, 
                                               # sheet = 'sampleinfo'
                                               sheet = 'sampleinfo.reorder'
                                               ),
                            stringsAsFactors = F)
    tumornormallist =  data.frame(readxl::read_excel(clinical.file, 
                                                     # sheet = 'tumor_normal'
                                                     sheet = 'tumor_normal.reorder'
                                                     ),
                                   stringsAsFactors = F)

    ## --------------------------------------------------------
    
    head(sampleinfo)
    # sampleinfo = unique(sampleinfo[,c('Sample','Type')])
    dim(sampleinfo) ## 61 11
    sampleinfo$Tumor_Sample_Barcode = sampleinfo$Sample
    dim(sampleinfo) ## 61 12
    tumornormallist$Tumor_Sample_Barcode = tumornormallist$Tumor
    mafs = NULL 
    
    ## import mafs
    for(my.sm in sort(tumornormallist$Tumor)) {
      
      print(my.sm)
      my.df = NULL 
      my.file = NULL 
      
      my.file = paste0('funcotator/',my.sm,
                       '.bwamem.mutect2.wpon.flt.flt.funcotator.fixed.maf')
      if(! file.exists(my.file)) { next }
      
      my.df = fread(my.file, stringsAsFactors = F, skip = 2860, sep = '\t')
      my.df = data.frame(my.df, stringsAsFactors = F)
      print(dim(my.df)) ## 366 222
      print(my.df[1:3,1:2])
      
      mafs = rbind(mafs,
                   my.df)
      
    }
    
    print(dim(mafs)) ## 23474   222
    table(mafs$Tumor_Sample_Barcode)
    length(table(mafs$Tumor_Sample_Barcode)) ## 52
    
    write.table(mafs,
                file = paste0(output,'.maf'), sep = '\t', row.names = F, quote = F)
    
    ## --------------------------------------------------------
    
    ## filter 
    dim(mafs) ## 23474   222
    
    if(FALSE) {
      # mafs = read.delim(paste0(output,'.fixed.maf'), stringsAsFactors = F,
      #                   na.strings = c('NA','na','','N/A','-','.'))
      class(mafs$gnomAD_exome_AF)
      mafs$gnomAD_exome_AF = gsub('[|]\\S+$','',mafs$gnomAD_exome_AF)
      mafs$gnomAD_exome_AF = gsub('[|]$','',mafs$gnomAD_exome_AF)
      mafs$gnomAD_exome_AF = as.numeric(as.character(mafs$gnomAD_exome_AF))
      
      mafs = mafs[is.na(mafs$gnomAD_exome_AF) | 
                    (!is.na(mafs$gnomAD_exome_AF) & (
                      mafs$gnomAD_exome_AF == '' | 
                        mafs$gnomAD_exome_AF < 0.0001
                    )),]
      dim(mafs) ## 272607    309
      
    }
      
     write.table(mafs,
                 file = paste0(output,'.fixed.gnomAD_exome_AF0.0001.maf'), sep = '\t', row.names = F, quote = F)
      
    
    
    ## --------------------------------------------------------
    
    mafs.0 = mafs
    basicStats(mafs.0$DP)
    # X..mafs.0.DP
    # nobs        1.023000e+04
    # NAs         0.000000e+00
    # Minimum     1.100000e+01
    # Maximum     3.328000e+03
    # 1. Quartile 4.100000e+01
    # 3. Quartile 2.640000e+02
    # Mean        1.899379e+02
    # Median      9.850000e+01
    # Sum         1.943065e+06
    # SE Mean     2.177730e+00
    # LCL Mean    1.856692e+02
    # UCL Mean    1.942067e+02
    # Variance    4.851586e+04
    # Stdev       2.202632e+02
    # Skewness    2.673704e+00
    # Kurtosis    1.518227e+01
    
    ## --------------------------------------------------------
    
    ## remove co-occuring sommut ... there are no known MLS mut hotspot! so likely to be artifact of unknown reason?!
    if(TRUE) {
      
    }
    
    
    ## --------------------------------------------------------
    
    for (my.dp in c(8,10,20)) {
      
      print(my.dp)
      
      ## ------------------------------------------------
      
      mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                       mafs.0$t_alt_count>=my.dp,]
      
      dim(mafs) ## 15113   309
      
      write.table(mafs,
                  file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'.maf'), 
                  sep = '\t', row.names = F, quote = F)
      
      ## ------------------------------------------------
      
      ## exclude dbsnp variants
      if(TRUE) {
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         mafs.0$dbSNP_RS=='',]
        dim(mafs) ## 4733  306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_exdbSNP.maf'), sep = '\t', row.names = F, quote = F)
        
        ## ------------------------------------------------
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         mafs.0$dbSNP_RS=='' &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f <=0.3,]
        dim(mafs) ## 831 306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_exdbSNP_altAF0.3max.maf'), sep = '\t', row.names = F, quote = F)
        
        ## ------------------------------------------------
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         mafs.0$dbSNP_RS=='' &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f >=0.1,]
        dim(mafs) ## 4502  306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_exdbSNP_altAF0.1min.maf'), sep = '\t', row.names = F, quote = F)
        
        ## ------------------------------------------------
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         mafs.0$dbSNP_RS=='' &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f <=0.3 &
                         mafs.0$tumor_f >=0.05,]
        dim(mafs) ## 735 306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_exdbSNP_altAF0.05min0.3max.maf'), sep = '\t', row.names = F, quote = F)
        
        ## ------------------------------------------------
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         mafs.0$dbSNP_RS=='' &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f <=0.3 &
                         mafs.0$tumor_f >=0.1,]
        dim(mafs) ## 600 306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_exdbSNP_altAF0.1min0.3max.maf'), sep = '\t', row.names = F, quote = F)
      }
      
      ## ------------------------------------------------
      
      ## NOT exclude dbsnp variants
      if(TRUE) {
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f <=0.3,]
        dim(mafs) ## 831 306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_altAF0.3max.maf'), sep = '\t', row.names = F, quote = F)
        
        ## ------------------------------------------------
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f >=0.1,]
        dim(mafs) ## 4502  306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_altAF0.1min.maf'), sep = '\t', row.names = F, quote = F)
        
        ## ------------------------------------------------
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f <=0.3 &
                         mafs.0$tumor_f >=0.05,]
        dim(mafs) ## 735 306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_altAF0.05min0.3max.maf'), sep = '\t', row.names = F, quote = F)
        
        ## ------------------------------------------------
        
        mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                         mafs.0$t_alt_count>=my.dp &
                         !is.na(mafs.0$tumor_f) &
                         mafs.0$tumor_f <=0.3 &
                         mafs.0$tumor_f >=0.1,]
        dim(mafs) ## 600 306
        
        write.table(mafs,
                    file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP',my.dp,'_altAF0.1min0.3max.maf'), sep = '\t', row.names = F, quote = F)
        
      }
      
    }
    
  }
  
  ## --------------------------------------------------------
  
  ## explore AF distribution ...
  if(TRUE) {
    
    my.dp = 8 ## 10 or 20, (9/16/2024) or 8.... 
    
    data.plot  = NULL 
    data.plot = read.delim(paste0(output,
                                  '.fixed.gnomAD_exome_AF0.0001_altDP',
                                  my.dp,
                                  '.maf'),
                           stringsAsFactors = F)
    
    dim(data.plot)
    
    p1 = NULL 
    p1 = ggplot(data.plot, aes(tumor_f)) +
      geom_density(aes(color = Tumor_Sample_Barcode)) +
      geom_vline(xintercept = c(0.3,0.4))
    
    p1
    
    p2 = NULL 
    p2 = ggplot(data.plot, aes(tumor_f)) +
      geom_histogram(aes(fill = Tumor_Sample_Barcode)) +
      # geom_vline(xintercept = c(0.3,0.4,0.7)) +
      facet_wrap(.~Tumor_Sample_Barcode, ncol = 6)
    
    p2
    
    ## decided to cap AF max as 0.3 ...
    
    ## one sample has a LOT of sommut: ....
    
  }
  
  ## --------------------------------------------------------
  
  ## rescue known driver gene mutations
  if(TRUE) {
    
    hnscc.known.genes = read.delim('HNSCC_HPVneg.known_sommut.txt', header = F)
    
    ## rescue TP53 and RB1
    if(TRUE) {
      
      out.file = paste0(output, '.fixed.gnomAD_exome_AF0.0001_altDP8_rescueKnownGenes.maf')
      
      in1.file = paste0(output, '.fixed.gnomAD_exome_AF0.0001_altDP8.maf')
      in2.file = paste0(output, '.fixed.gnomAD_exome_AF0.0001_altDP8_exdbSNP.maf')
      
      in1 = read.delim(in1.file)
      in2 = read.delim(in2.file)
      
      # has mut in known genes and also have dbSNP id
      in1 = in1[in1$Hugo_Symbol %in% hnscc.known.genes$V1 &
                  !(is.na(in1$dbSNP_ID) | in1$dbSNP_ID ==''),,drop=F]
      dim(in1) ## 17 222
      
      dim(in2) ## 3062  222
      in2 = rbind(in1, in2)
      dim(in2) ##  3079  222
      
      write.table(in2, 
                  file = out.file,
                  sep = '\t', row.names = F, col.names = T, quote = F)
      
    }
   
    ## --------------------------------------------------------
    
    
    ## rescue TP53 and RB1, also exclude genes that are mut too frequently...
    if(FALSE) {
      
      out.file = paste0(output, '.fixed.gnomAD_exome_AF0.0001_altDP10_exdbSNP_rescueKnownGenes_altAF0.3max_geneMutFrac0.1max.maf')

      in1.file =paste0(output, '.fixed.gnomAD_exome_AF0.0001_altDP10_exdbSNP_rescueKnownGenes_altAF0.3max.maf')
      
      in1 = read.delim(in1.file)

      x = NULL 
      x = in1[in1$Variant_Classification %in% var.funcs,
              c('Hugo_Symbol','Tumor_Sample_Barcode')]
      dim(x) ## 7453    2
      x = unique(x)
      dim(x) ## 7007    2
      x = data.frame(table(x$Hugo_Symbol))
      x$Var1 = as.character(x$Var1)
      x = x[rev(order(x$Freq)),]
      x$Total = length(unique(in1$Tumor_Sample_Barcode))
      x$Frac = x$Freq / x$Total
      
      dim(x) ## 5056    4
      dim(x[x$Frac >= 0.1,])
      ## 22 
      # Var1 Freq Total      Frac
      # 4577     TTN   14    56 0.2500000
      # 4946  ZNF493   11    56 0.1964286
      # 2651   MUC16   10    56 0.1785714
      # 2220   KMT2D   10    56 0.1785714
      # 3730    RYR1    9    56 0.1607143
      # 3666  RNF213    7    56 0.1250000
      # 3139    PCLO    7    56 0.1250000
      # 1931   HSPG2    7    56 0.1250000
      # 1213   DNAH8    7    56 0.1250000
      # 1032   CSMD3    7    56 0.1250000
      # 12    ABCA13    7    56 0.1250000
      # 3493   PTPRF    6    56 0.1071429
      # 3239   PKHD1    6    56 0.1071429
      # 2839  NOTCH4    6    56 0.1071429
      # 2652   MUC17    6    56 0.1071429
      # 2618  MROH2B    6    56 0.1071429
      # 1887   HMCN1    6    56 0.1071429
      # 1485  FAM83H    6    56 0.1071429
      # 1451 FAM135B    6    56 0.1071429
      # 935  COL24A1    6    56 0.1071429
      # 672  CCDC168    6    56 0.1071429
      # 140   ADGRB1    6    56 0.1071429
      
      
      print(x[x$Frac>=0.1 & x$Var1 %in% hnscc.known.genes$V1,])
      # Var1 Freq Total      Frac
      # 2220  KMT2D   10    56 0.1785714
      # 2839 NOTCH4    6    56 0.1071429
      ## distribution here for the two genes is a bit weird ....
      ## https://aacrjournals.org/clincancerres/article/25/19/5961/124617/The-Genomic-Landscape-of-Merkel-Cell-Carcinoma-and
      
      x[x$Var1 %in% c('TP53','RB1'),]
      # Var1 Freq Total       Frac
      # 4459 TP53    1    56 0.01785714
      # 3560  RB1    1    56 0.01785714
      
      ## I decided to exclude all genes with mut prevelance of 10% or higher 
      x = x[x$Frac >=0.1,,drop=F]
      dim(x) # 22
      
      dim(in1) ## 23469   309
      in1 = in1[!in1$Hugo_Symbol %in% x$Var1,]
      dim(in1) ## 23033   309
      
      write.table(in1, 
                  file = out.file,
                  sep = '\t', row.names = F, col.names = T, quote = F)
    }
    
  }
}

## --------------------------------------------------------

## data files 
if(TRUE) {
  
  hnscc.known.genes.file = 'HNSCC_HPVneg.known_sommut.txt'
  hnscc.known.genes = read.delim(hnscc.known.genes.file, header = F)
  
  ## --------------------------------------------------------
  
  clinical.file = '../../../sampleinfo/PDO_sm61.exome.sampleinfo.20250102.xlsx'
  
  sampleinfo = data.frame(readxl::read_excel(clinical.file, 
                                             # sheet = 'sampleinfo'
                                             sheet = 'sampleinfo.reorder'
  ),
  stringsAsFactors = F)
  tumornormallist =  data.frame(readxl::read_excel(clinical.file, 
                                                   # sheet = 'tumor_normal'
                                                   sheet = 'tumor_normal.reorder'
  ),
  stringsAsFactors = F)
  
  clinical = tumornormallist
  
  ## --------------------------------------------------------
  
  maf.file = paste0(output, '.fixed.gnomAD_exome_AF0.0001.maf')
  # maf.file.2 = paste0(output, '.fixed.gnomAD_exome_AF0.0001_altDP8.maf')
  # maf.file.3 = paste0(output, '.fixed.gnomAD_exome_AF0.0001_altDP8_rescueKnownGenes.maf')

}

## --------------------------------------------------------

in.file = NULL 
in.file = maf.file
# in.file = maf.file.2
# in.file = maf.file.3

## --------------------------------------------------------

## explore  TMB as a result of each maf file...
## compute TMB, add known mutations
if(TRUE) {
 
  
  x = NULL 
  x = read.delim(in.file, stringsAsFactors = F)
  
  dim(x) ##  23474   222
  x = x[!x$Tumor_Sample_Barcode %in% sample.ex,]
  head(x$dbSNP_ID)
  
  ## how much % of variants have snpid??
  print(nrow(x[!is.na(x$dbSNP_ID) & !x$dbSNP_ID %in% c(''),,drop=F]) / 
          nrow(x))
  ## 0.2262503
  
  dim(x) # 23474   222
  # x = x[is.na(x$dbSNP_ID) | (x$dbSNP_ID=='') | 
  #         (!is.na(x$dbSNP_ID) & x$dbSNP_ID!='' & x$Hugo_Symbol %in% c('TP53')),]
  # dim(x)
  
  print(length(unique(x$Tumor_Sample_Barcode))) ## 52
  
  ## --------------------------------------------------------
  
  exome.tmb = NULL 
  exome.tmb = data.frame(Sample.exome = sort(unique(x$Tumor_Sample_Barcode)),
                          stringsAsFactors = F)
  
  ## --------------------------------------------------------
  
  table(x$Variant_Classification[x$Tumor_Sample_Barcode %in% c('TIIL_0165')])
  # 3'UTR             5'Flank               5'UTR COULD_NOT_DETERMINE 
  #                18                   2                  12                   2 
  #   Frame_Shift_Del        In_Frame_Del              Intron   Missense_Mutation 
  #                 2                   3                 191                  66 
  # Nonsense_Mutation                 RNA              Silent         Splice_Site 
  #                11                   7                  33                  19 
  
     
     
  ## coding variants only
  x = x[x$Variant_Classification %in% var.funcs,]
  dim(x) #  5653  222
  
  print(length(unique(x$Tumor_Sample_Barcode))) ## 52
  ## 1 sample filtered out??
  
  print(exome.tmb[!exome.tmb$Sample.exom %in% x$Tumor_Sample_Barcode,,drop=F])
## all in!
  
  write.table(x,
              file = paste0(in.file,'.coding.maf'),
              sep = '\t', quote = F, row.names = F, col.names = T)
  
  
  ## --------------------------------------------------------
  
  y = NULL 
  y = data.frame(table(x$Tumor_Sample_Barcode))
  colnames(y) = c('Sample.exome','TMB')
  
  y
  
  exome.tmb = merge(exome.tmb, y, by = 'Sample.exome',all.x = T)
  exome.tmb[is.na(exome.tmb)] = 0
  dim(exome.tmb) ## 25 samples 
  
  p1 = ggplot(exome.tmb, aes(TMB)) +
    geom_density()
  
  p1
  
  ## 600 might be the high TMB cutoff?!
  median(exome.tmb$TMB)
  ## 90
  
  exome.tmb[exome.tmb$TMB >=600,]
  ### none!
  
  exome.tmb$TMB.perMB = exome.tmb$TMB / 50
  exome.tmb[exome.tmb$TMB.perMB >=10,]
  # Sample.exome TMB TMB.perMB
  # 10    TIIL_0177 556     11.12
  
  x[x$Hugo_Symbol %in% c('TP53'),
    c('Tumor_Sample_Barcode','Hugo_Symbol','tumor_f','Protein_Change','cDNA_Change')]
  # Tumor_Sample_Barcode Hugo_Symbol tumor_f Protein_Change cDNA_Change
  # 283             TIIL_0165        TP53   0.220        p.R248W    c.742C>T
  # 647             TIIL_0167        TP53   0.989        p.R248W    c.742C>T
  # 1011            TIIL_0168        TP53   0.993        p.R248W    c.742C>T
  # 1321            TIIL_0169        TP53   0.226        p.R175H    c.524G>A
  # 3399            TIIL_0177        TP53   0.380        p.R337L   c.1010G>T
  # 3400            TIIL_0177        TP53   0.329        p.T125T    c.375G>T
  # 4443            TIIL_0180        TP53   0.628                  c.e6-1G>A
  # 4952            TIIL_0182        TP53   0.994                  c.e6-1G>A
  # 5422            TIIL_0183        TP53   0.993                  c.e6-1G>A
  # 6012            TIIL_0184        TP53   0.985                  c.e6-1G>A
  # 6532            TIIL_0185        TP53   0.992                  c.e6-1G>A
  # 7426            TIIL_0186        TP53   0.088        p.R248W    c.742C>T
  # 7427            TIIL_0186        TP53   0.073        p.V203L    c.607G>T
  # 7428            TIIL_0186        TP53   0.066        p.A161T    c.481G>A
  # 8842            TIIL_0190        TP53   0.990        p.R248W    c.742C>T
  # 9628            TIIL_0193        TP53   0.981                  c.e6-1G>A
  ## ..... more 
  
  x[x$Hugo_Symbol %in% hnscc.known.genes$V1,
    c('Tumor_Sample_Barcode','Hugo_Symbol','tumor_f','Protein_Change','cDNA_Change')]
  
  exome.tmb
  # Sample.exome TMB TMB.perMB
  # 1     TIIL_0165 101      2.02
  # 2     TIIL_0167  96      1.92
  # 3     TIIL_0168 103      2.06
  # 4     TIIL_0169  50      1.00
  # 5     TIIL_0171   5      0.10
  # 6     TIIL_0172  15      0.30
  # 7     TIIL_0173  12      0.24
  # 8     TIIL_0174  23      0.46
  # 9     TIIL_0176  23      0.46
  # 10    TIIL_0177 556     11.12
  # 11    TIIL_0179  24      0.48
  # 12    TIIL_0180  82      1.64
  # 13    TIIL_0182  94      1.88
  # 14    TIIL_0183  86      1.72
  # 15    TIIL_0184 108      2.16
  # 16    TIIL_0185  88      1.76
  # 17    TIIL_0186 333      6.66
  # 18    TIIL_0188  31      0.62
  # 19    TIIL_0189  52      1.04
  # 20    TIIL_0190 104      2.08
  # 21    TIIL_0191  19      0.38
  # 22    TIIL_0192  15      0.30
  # 23    TIIL_0193  81      1.62
  # 24    TIIL_0194  22      0.44
  # 25    TIIL_0195  24      0.48
  # 26    TIIL_0260 105      2.10
  # 27    TIIL_0261  73      1.46
  # 28    TIIL_0262 121      2.42
  # 29    TIIL_0264 133      2.66
  # 30    TIIL_0265 139      2.78
  # 31    TIIL_0266 146      2.92
  # 32    TIIL_0267 132      2.64
  # 33    TIIL_0268 134      2.68
  # 34    TIIL_0269 143      2.86
  # 35    TIIL_0270 139      2.78
  # 36    TIIL_0271  22      0.44
  # 37    TIIL_0273  80      1.60
  # 38    TIIL_0274  83      1.66
  # 39    TIIL_0275  87      1.74
  # 40    TIIL_0276  87      1.74
  # 41    TIIL_0277  83      1.66
  # 42    TIIL_0278  79      1.58
  # 43    TIIL_0279  79      1.58
  # 44    TIIL_0280 142      2.84
  # 45    TIIL_0282  92      1.84
  # 46    TIIL_0283 233      4.66
  # 47    TIIL_0284 195      3.90
  # 48    TIIL_0285 131      2.62
  # 49    TIIL_0287 299      5.98
  # 50    TIIL_0288 184      3.68
  # 51    TIIL_0289 185      3.70
  # 52    TIIL_0290 180      3.60
  
  ## --------------------------------------------------------
  
  z = NULL 
  z = x[x$Hugo_Symbol %in% hnscc.known.genes$V1,,drop=F]
  dim(z) ## 76 222
  z
  
  write.table(z,
              file = paste0(in.file,'.coding.known.maf'),
              sep = '\t', quote = F, row.names = F, col.names = T)
  
  ## --------------------------------------------------------
  
  z = z[,c('Tumor_Sample_Barcode','Hugo_Symbol','Protein_Change','cDNA_Change','Variant_Classification'),
        drop=F]
  
  row.select = NULL 
  row.select = which(z$Protein_Change=='')
  z$Protein_Change[row.select] = z$Variant_Classification[row.select]
  
  z = reshape2::dcast(Tumor_Sample_Barcode ~ Hugo_Symbol,
                      value.var = 'Protein_Change',
                      data = z,
                      fun.aggregate = function(x) { paste(x, collapse = ',')})
  
  row.names(exome.tmb) = exome.tmb$Sample.exome
  row.names(z) = z$Tumor_Sample_Barcode
  
  dim(exome.tmb)
  exome.tmb = merge(exome.tmb, z[,-1], by = 'row.names', all.x = T)
  dim(exome.tmb) ## 25  8
  
  row.names(exome.tmb) = exome.tmb$Row.names
  exome.tmb = exome.tmb[,-1]
  
  ## --------------------------------------------------------
  
  exome.tmb[exome.tmb==''] = NA
  
  write.table(exome.tmb,
              file = paste0(in.file,'.coding.TMB.txt'),
              sep = '\t', quote = F, row.names = F, col.names = T)
  
  write.csv(exome.tmb,
              file = paste0('HNSCC.exome.sm',nrow(exome.tmb),'.tmb.csv'))
  
}

## --------------------------------------------------------

## import data 
if(FALSE) {
  
  # # hnscc.known.genes = read.delim(hnscc.known.genes.file, header = F)
  # 
  # 
  # clinical = data.frame(readxl::read_excel(clinical.file, sheet = 'merged',
  #                                          na = c('NA','','N/A','na','n/a')),
  #                         stringsAsFactors = F)
  # rnaseq.viral =  data.frame(readxl::read_excel(clinical.file, sheet = 'giahsk.kraken2.c0.1',
  #                                               na = c('NA','','N/A','na','n/a')),
  #                            stringsAsFactors = F)
  # clinical$Tumor_Sample_Barcode = clinical$Sample.exome
  # 
  # ## keep only entries with exome data avail
  # dim(clinical) ## 91 114
  # clinical = clinical[!is.na(clinical$Sample.exome),]
  # dim(clinical) ## 56 114
  
  maf = read.maf(in.file, clinicalData = clinical)
  if(length(sample.ex) > 0) {
    maf = filterMaf(maf, tsb = sample.ex)
    
  }
  
}

## --------------------------------------------------------

## draw plots outside of maftools 
if(FALSE) {
  
  rm(x,y,z)
  
  x = read.delim(maf.file,stringsAsFactors = F)
  dim(x) # 762 309
  
  table(x$Variant_Classification)
  
  ## all sommut 
  if(TRUE) {
    
    y = data.frame(table(x[,c('Tumor_Sample_Barcode','Variant_Classification')]))
    print(sum(y$Freq)) ## 34428
    
    library(RColorBrewer)
    set.seed(10)
    n <- length(unique(x$Variant_Classification))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    n.colors = sample(col_vector, n)
    pie(rep(1,n), col=n.colors)
    
    
    plot.order = y
    plot.order = reshape2::dcast(Tumor_Sample_Barcode ~ Variant_Classification,
                                 data = y,
                                 value.var = 'Freq',
                                 fun.aggregate = sum)
    print(plot.order)
    
    write.csv(plot.order,
              file = paste0(maf.file,'.var',nrow(x),'.csv'))
    
    row.names(plot.order) = plot.order$Tumor_Sample_Barcode
    plot.order = data.frame(total = apply(plot.order[,-1], 1, sum))
    plot.order = plot.order[rev(order(plot.order$total)),,drop=F]
    
    y$Tumor_Sample_Barcode = factor(y$Tumor_Sample_Barcode,
                                    levels = as.character(row.names(plot.order)))
    
    
    p1 = ggplot(y, aes(Variant_Classification, Freq)) +
      geom_bar(aes(fill = Variant_Classification), width = 0.6, 
               position = 'stack', stat = 'identity', color = '#000000') +
      scale_fill_manual(values = n.colors) +
      theme_pubr(x.text.angle = 90) +
      facet_grid(Tumor_Sample_Barcode~.)
    
    p1
    
    pdf(file = paste0(maf.file,'.var',nrow(x),'.pdf'),
        width = 12, height = 12)
    print(p1)
    dev.off()
    
  }
  
  ## coding variants only
  if(TRUE) {
    
    y = data.frame(table(x[x$Variant_Classification %in% 
                             c(var.funcs,'Silent'),
                           c('Tumor_Sample_Barcode','Variant_Classification')]))
    print(sum(y$Freq)) ## 181
    
    library(RColorBrewer)
    set.seed(10)
    n <- length(unique(x$Variant_Classification))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    n.colors = sample(col_vector, n)
    pie(rep(1,n), col=n.colors)
    
    
    plot.order = y
    plot.order = reshape2::dcast(Tumor_Sample_Barcode ~ Variant_Classification,
                                 data = y,
                                 value.var = 'Freq',
                                 fun.aggregate = sum)
    print(plot.order)
    
    write.csv(plot.order,
              file = paste0(maf.file,'.var',nrow(x[x$Variant_Classification %in% 
                                                     c(var.funcs,'Silent'),,drop=F]),'.csv'))
    
    row.names(plot.order) = plot.order$Tumor_Sample_Barcode
    plot.order = data.frame(total = apply(plot.order[,-1], 1, sum))
    plot.order = plot.order[rev(order(plot.order$total)),,drop=F]
    
    y$Tumor_Sample_Barcode = factor(y$Tumor_Sample_Barcode,
                                    levels = as.character(row.names(plot.order)))
    
    
    p1 = ggplot(y, aes(Variant_Classification, Freq)) +
      geom_bar(aes(fill = Variant_Classification), width = 0.6, 
               position = 'stack', stat = 'identity', color = '#000000') +
      scale_fill_manual(values = n.colors) +
      theme_pubr(x.text.angle = 90) +
      facet_grid(Tumor_Sample_Barcode~.)
    
    p1
    
    pdf(file = paste0(maf.file,'.var',nrow(x[x$Variant_Classification %in% 
                                               c(var.funcs,'Silent'),,drop=F]),'.pdf'),
        width = 12, height = 12)
    print(p1)
    dev.off()
    
  }
 
}

## --------------------------------------------------------

## check existing mutations (somatic? germline?) from literature 
if(FALSE) {
  
  my.genes = c('PIK3CA','PTEN','KMT2D','EP300','MTOR','SDHA','TP53','AKT1','KMT2A','ARID1B','PDGFRA','FLT1','FGFR1','FANCA','RECQL4','DNMT3A','ATM')
  
  my.df = read.delim('MLS.exome.fixed.gnomAD_exome_AF0.0001_altDP5.maf',stringsAsFactors = F)
  table(my.df$Variant_Classification)
  
  
  my.df = my.df[my.df$Hugo_Symbol %in% my.genes &
                  my.df$Variant_Classification %in% 
                  c('DE_NOVO_START_IN_FRAME','DE_NOVO_START_OUT_FRAME',
                    'Frame_Shift_Del','Frame_Shift_Ins',
                    'In_Frame_Del','In_Frame_Ins',
                    'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation',
                    'Splice_Site',
                    'START_CODON_SNP'),]
  dim(my.df)
  
  write.csv(my.df,
            file = paste0('MLS.exome.fixed.gnomAD_exome_AF0.0001_altDP5.maf',
                          '.known_MLS_genes.csv'))
  
}

## --------------------------------------------------------

## parse maf
if(FALSE) {
  
  maf.sampleSummary = getSampleSummary(maf)
  median(maf.sampleSummary$total)
  print(maf.sampleSummary)
  
  basicStats(maf.sampleSummary$total / 35)
  
  clinical[!clinical$Tumor_Sample_Barcode %in% maf.sampleSummary$Tumor_Sample_Barcode,] ## all samples in!!
  
  write.csv(maf.sampleSummary,
            file = paste0(in.file,'.getSampleSummary.csv'))
  
  pdf(file = paste0(in.file,'.landscape.pdf'), width = 10, height = 8)
  
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', 
                 dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  pdf(file = paste0(in.file,'.oncoprint.pdf'), width = 10, height = 6)
  
  ## oncoplot for top ten mutated genes.
  oncoplot(maf = maf, top = 20, removeNonMutated = F,
           clinicalFeatures = c('PID'),
           showTumorSampleBarcodes = T,
           sampleOrder = clinical$Tumor_Sample_Barcode)
  dev.off()
  
  
  maf.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
  #plot titv summary
  plotTiTv(res = maf.titv)
  
  plotVaf(maf = maf, vafCol = 'tumor_f')
  
  #exclusive/co-occurance event analysis on top 10 mutated genes. 
  somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
  
  ## Detecting cancer driver genes based on positional clustering
  maf.sig = oncodrive(maf = maf, AACol = 'Protein_Change', 
                      minMut = 5, pvalMethod = 'zscore')
  head(maf.sig)
  
  plotOncodrive(res = maf.sig, fdrCutOff = 0.15, useFraction = TRUE)
  rm(maf.sig)

  lollipopPlot(maf = maf, 
               gene = 'TP53', labelPos = 'all',
               AACol = 'Protein_Change', showMutationRate = TRUE)
  
  ## https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01068-0
  # We nominated PTPRJ, mutated and homozygously deleted in 3.8% (9/240) and 0.8% (2/240) of samples, respectively, as a probable tumor suppressor gene, and FER and SKP2, amplified in 3.8% and 11.7% of samples, respectively, as probable oncogenes. We further identified a long tail of infrequent pathogenic alterations, involving genes such as CIC and LZTR1. Pathogenic germline mutations were observed on MITF, PTEN, ATM, and PRKN. We found BRAF V600E mutations in acral melanomas with fewer structural variations, suggesting that they are distinct and related to cutaneous melanomas. Amplifications of PAK1 and GAB2 were more commonly observed in acral melanomas, whereas SF3B1 R625 codon mutations were unique to mucosal melanomas (12.9%). Amplifications at 11q13-14 were frequently accompanied by fusion to a region on chromosome 6q12, revealing a recurrent novel structural rearrangement whose role remains to be elucidated.
  
  
  ## Adding and summarizing pfam domains
  pfamDomains(maf = maf, AACol = 'Protein_Change', top = 10)
  
  ## total TMB 
  maf.samplesummary = getSampleSummary(maf)
  clinical$Tumor_Sample_Barcode = as.character(clinical$Tumor_Sample_Barcode)
  maf.samplesummary$Tumor_Sample_Barcode = as.character(maf.samplesummary$Tumor_Sample_Barcode)
  
  clinical = merge(clinical, maf.samplesummary, by = 'Tumor_Sample_Barcode')
  
  library(fBasics)
  basicStats(clinical$total)
  

  aggregate(total ~ PID ,
            data = clinical,
            FUN = mean)
  
  aggregate(total ~ PID ,
            data = clinical,
            FUN = std.error)
  
  write.csv(clinical,
            file = paste0(clinical.file,'.wSeqfiles.wTMB.csv'))
  
  p1 = ggplot(clinical, aes(PID, total)) +
    geom_boxplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 5, width = 0.1, height = 0, fill = '#999999') +
    # scale_fill_manual(values = plot.colors) + 
    # stat_compare_means(method = 'wilcox.test') +
    ggtitle('Total TMB (NSSMs)') +
    theme_classic() +
    stat_compare_means()
  
  p1
  
  ## log transform ...
  min(clinical$total) ## 2
  max(clinical$total) ## 337
  dim(clinical) ## 40 samples 
  
  p2.1 = ggplot(clinical, aes(PID, total)) +
    geom_boxplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 5, width = 0.1, height = 0, fill = '#999999') +
    # scale_fill_manual(values = plot.colors) + 
    stat_compare_means(method = 't.test') +
    ggtitle('Total TMB (Non-Synonymous Somatic Mutations (SNVs/indels), NSSMs)') +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = 'l')  +
    theme_pubr()
  
  p2.1
  
  p3 = ggplot(clinical, aes(total)) +
    geom_density()
  
  p3
  
  p3.1 = ggplot(clinical, aes(log10(total))) +
    geom_density()
  
  p3.1
  
  
}

## --------------------------------------------------------

## run mut signature 
if(FALSE) {

  cosmic.sig.db ="SBS" ## SBS (65) legacy (30)
  
  ## --------------------------------------------------
  
  print(maf.file)

  ## --------------------------------------------------
  
  ## mutsig detection method 1
  if(TRUE) {
    
    ## only works in R 3.6.0 .... does not work in R 4.0.3 
    
    library(deconstructSigs)  
    library(BSgenome.Hsapiens.UCSC.hg38)
    
    sample.mut.ref = read.delim(maf.file, stringsAsFactors = F)

    head(sample.mut.ref)
    
    # Convert to deconstructSigs input
    sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, 
                                    sample.id = "Tumor_Sample_Barcode", 
                                    chr = "Chromosome", 
                                    pos = "Start_Position", 
                                    ref = "Reference_Allele", 
                                    alt = "Tumor_Seq_Allele2",
                                    bsg = BSgenome.Hsapiens.UCSC.hg38)
    #> Creating a generic function for 'nchar' from package 'base' in package 'S4Vectors'
    
    dim(sigs.input)
    sigs.input[1:4,1:3]
    
    ## --------------------------------------------------
    
    for(analysis in 1:4) {
      
      out.dir = NULL 
      sig.ref = NULL 
      tri.method = NULL 
      
      if(analysis == 1) {
        sig.ref = signatures.nature2013
        tri.method = 'default'
      } else if (analysis == 2) {
        sig.ref = signatures.nature2013
        tri.method = 'genome'
        
      } else if(analysis == 3) {
        sig.ref = signatures.cosmic
        tri.method = 'default'
      } else if (analysis == 4) {
        sig.ref = signatures.cosmic
        tri.method = 'genome'
      }
      
      ## --------------------------------------------------
      
      out.dir = file.path('deconstructSigs',
                          paste0('var',nrow(sample.mut.ref)),
                          paste0('sig',nrow(sig.ref),
                                 '_tri',tri.method))
      
      if( !dir.exists(out.dir)  ){ dir.create(out.dir, recursive = T) }
      
      ## --------------------------------------------------
      
      sigs.output = list()
      sigs.weights = NULL 
      
      for(i in 1:nrow(sigs.input)) {
        
        my.sm = row.names(sigs.input)[i]
        print(my.sm)
        my.df = NULL 
        
        # Determine the signatures contributing to the two example samples
        my.df = whichSignatures(tumor.ref = sigs.input, 
                                signatures.ref = sig.ref, 
                                sample.id = my.sm, 
                                contexts.needed = TRUE,
                                tri.counts.method = tri.method)
        
        sigs.output[[i]] = my.df
        sigs.weights = rbind(sigs.weights,
                             my.df$weights)
        
        ## --------------------------------------------------
        
        my.out = NULL 
        my.out = file.path(out.dir,
                           paste0(maf.file,
                                  '.deconstructSigs.',my.sm))
        
        saveRDS(my.df,
                file = paste0(my.out, '.rds'))
        
        write.csv(my.df$weights,
                  file = paste0(my.out, '.csv'))
        
        # Plot output
        pdf(file = paste0(my.out,'.bar.pdf'), height = 5, width = 10)
        print(plotSignatures(my.df, sub = my.sm))
        dev.off()
        
        
        pdf(file = paste0(my.out,'.pie.pdf'), height = 5, width = 10)
        print(makePie(my.df, sub = my.sm))
        dev.off()
        
      }
      
      ## --------------------------------------------------
      
      my.out = file.path(out.dir,
                         paste0(maf.file,
                                '.deconstructSigs.sm',nrow(sigs.input)))
      
      dim(sigs.weights)
      write.csv(sigs.weights,
                file = paste0(my.out,'.csv'))
      
      ## --------------------------------------------------
      
      data.plot = NULL 
      data.plot = data.frame(Sample = row.names(sigs.weights),
                             sigs.weights,
                             stringsAsFactors = F)
      data.plot = reshape2::melt(data.plot)
      colnames(data.plot)[2:3] = c('sig.ref','weights')
      
      p1 = NULL 
      p1 = ggplot(data.plot, aes(sig.ref, weights)) +
        geom_bar(stat = 'identity') +
        facet_grid(Sample~.) +
        theme_pubr(x.text.angle = 90)
      
      p1
      
      pdf(paste0(my.out,'.bar.pdf'), height = 20, width = 20)
      print(p1)
      dev.off()
      
      
      
    }
    
    
  }
  
  ## --------------------------------------------------
  
  ## mutsig detection method 2
  if(TRUE) {
    
    library(sigminer)
    
    ## import data
    if(TRUE) {
      
      # maf = read.delim(maf.file, stringsAsFactors = F)
      # dim(maf) ## 2621  307
      
      cosmic.sig.db ="SBS"  ## SBS DBS ID 
      
      # maf =  read_maf(maf = maf.file)
      
      # dim(maf)
      
      # maf.input = maf[,c('Tumor_Sample_Barcode','Chromosome',
      #                    'Start_Position','End_Position',
      #                    'Reference_Allele','Tumor_Seq_Allele2')]
      
      rm(x)
      x = read.delim(maf.file,stringsAsFactors = F)
      dim(x) # 762 309
      
      sample.mut.ref = x
      head(sample.mut.ref)
      
      # Convert to deconstructSigs input
      sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, 
                                      sample.id = "Tumor_Sample_Barcode", 
                                      chr = "Chromosome", 
                                      pos = "Start_Position", 
                                      ref = "Reference_Allele", 
                                      alt = "Tumor_Seq_Allele2",
                                      bsg = BSgenome.Hsapiens.UCSC.hg38)
      
      
      
      # Determine the signatures contributing an already normalized sample
      test = whichSignatures(tumor.ref = randomly.generated.tumors, 
                             signatures.ref = signatures.nature2013, 
                             sample.id = 2)
      
    }
    
    ## --------------------------------------------------
    
    ## de novo detection of signature
    if(TRUE) {
      
      
      
      mt_tally <- sig_tally(
        maf,
        ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
        useSyn = TRUE
      )
      
      mt_tally$nmf_matrix[1:5, 1:5]
      
      samples=sort(unique(row.names( mt_tally$nmf_matrix)))
      
      if(TRUE) {
        
        for(my.sm in samples) {
          print(my.sm)
          
          p1 = NULL 
          p1 = show_catalogue(t(mt_tally$nmf_matrix), style = "cosmic", x_label_angle = 90,
                              samples_name = my.sm,
                              samples = my.sm)
          pdf(file = paste0(maf.file,'.',my.sm,'.sigminer.raw.pdf'),
              height = 3, width = 10)
          print(p1)
          dev.off()
          
        }
      }
      
      
      
      # ## --------------------------------------------------
      # 
      # mt_tally_all <- sig_tally(
      #   maf,
      #   ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
      #   useSyn = TRUE,
      #   mode = "ALL",
      #   add_trans_bias = TRUE
      # )
      # 
      # str(mt_tally_all, max.level = 1)
      
    }
    
    ## --------------------------------------------------
    
    ## extract signature 
    if(TRUE) {
      
      ## in this ESC WGS cohort, any NMF after 4 returned error ...
      # Compute NMF rank= 2  ... + measures ... OK
      # Compute NMF rank= 3  ... + measures ... OK
      # Compute NMF rank= 4  ... + measures ... OK
      # Compute NMF rank= 5  ... + measures ... ERROR
      # Compute NMF rank= 6  ... + measures ... ERROR
      # Compute NMF rank= 7  ... + measures ... ERROR
      # Compute NMF rank= 8  ... + measures ... ERROR
      # Compute NMF rank= 9  ... + measures ... ERROR
      # Compute NMF rank= 10  ... Timing stopped at: 0 0 0 (I manually killed the job)
      
      mt_est <- sig_estimate(mt_tally$nmf_matrix,
                             # range = 2:10,
                             range = 2:4,
                             nrun = 10, # increase this value if you wana a more stable estimation
                             use_random = FALSE, # if TRUE, add results from randomized input
                             cores = 4,
                             # pConstant = 1e-13,
                             verbose = TRUE
      )
      
      
      ## You can also select the measures to show
      ## by 'what' option
      show_sig_number_survey2(mt_est$survey)
      
      show_sig_number_survey(mt_est$survey)
      
      ## --------------------------------------------------
      
      mt_sig <- sig_extract(mt_tally$nmf_matrix,
                            # n_sig = 6,
                            n_sig = 4,
                            nrun = 30,
                            cores = 4
                            # pConstant = 1e-13
      )
      
      sim_v3 <- get_sig_similarity(mt_sig, sig_db = cosmic.sig.db)
      
      pheatmap::pheatmap(sim_v3$similarity)
      
      ## --------------------------------------------------
      
      sig_signature(mt_sig)[1:4, ]
      
      sig_exposure(mt_sig)[, 1:5]
      
      get_sig_exposure(mt_sig)
      
    }
    
    ## --------------------------------------------------
    
    ## make plots
    if(TRUE) {
      
      my.sm = samples[1:3]
      print(my.sm)
      show_catalogue(t(mt_tally$nmf_matrix), style = "cosmic", x_label_angle = 90,
                     samples_name = samples,
                     samples = samples)
      
      ## --------------------------------------------------
      
      show_catalogue(t(mt_tally$nmf_matrix), style = "cosmic", x_label_angle = 90)
      
      # show_catalogue(t(mt_tally$nmf_matrix), style = "cosmic", samples = "TCGA-E9-A22E-01A-11D-A159-09", x_label_angle = 90)
      
      
      ## --------------------------------------------------
      
      dim(mt_sig$Exposure)
      
      mt_sig.exp = mt_sig$Exposure
      mt_sig.exp.norm = mt_sig$Exposure.norm
      
      ## --------------------------------------------------
      
      show_sig_exposure(mt_sig, style = "cosmic",
                        hide_samps = F)
      
      # ## exome 
      # show_sig_exposure(mt_sig, style = "cosmic",
      #                   sig_names = c('S1:DDMR_clock_chemo',
      #                                 'S2:DDMR_clock_chemo',
      #                                 'S3:clock',
      #                                 'S4:DDMR_clock_chemo',
      #                                 'S5:clock',
      #                                 'S6:clock'),
      #                   hide_samps = F)
      # 
      # ## genome 
      # show_sig_exposure(mt_sig, style = "cosmic",
      #                   sig_names = c('S1:clock',
      #                                 'S2:APOBEC',
      #                                 'S3:clock',
      #                                 'S4:DDMR_clock_chemo',
      #                                 'S5:HRD_clock',
      #                                 'S6:Chemo_HRD_clock'),
      #                   hide_samps = F)
      
      # ## exome_genome 
      # show_sig_exposure(mt_sig, style = "cosmic",
      #                   sig_names = c('S1:clock',
      #                                 'S2:APOBEC',
      #                                 'S3:clock',
      #                                 'S4:DDMR_clock_chemo',
      #                                 'S5:HRD_clock',
      #                                 'S6:Chemo_HRD_clock'),
      #                   hide_samps = F)
      # 
      # show_sig_exposure(sig_w, style = "cosmic")
      
      # grp <- get_groups(sig_w)
      # 
      # grp_label <- grp$group
      # names(grp_label) <- grp$sample
      # 
      # show_sig_exposure(sig_w, style = "cosmic", groups = grp_label)
      # 
      ## --------------------------------------------------
      
      show_sig_consensusmap(mt_sig)
      
      
    }
    
  }
  
}

## --------------------------------------------------------

## compare complex vs isolated karyote groups
if(FALSE) {
  
  print(maf.file)
  
  maf = read.delim(maf.file, stringsAsFactors = F)
  
  table(maf$Variant_Classification)
  
  
  maf = maf[maf$Variant_Classification %in% 
                  c('DE_NOVO_START_IN_FRAME','DE_NOVO_START_OUT_FRAME',
                    'Frame_Shift_Del','Frame_Shift_Ins',
                    'In_Frame_Del','In_Frame_Ins',
                    'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation',
                    'Splice_Site',
                    'START_CODON_SNP'),]
  
  dim(maf) ## 162  309
  
  ## --------------------------------------------------------
  
  group1 = paste0('S',sort(unique(clinical$`CBL ID`[clinical$Karyotype...5=='complex'])))
  group2 = paste0('S',sort(unique(clinical$`CBL ID`[clinical$Karyotype...5=='Isolated'])))
  
  print(group1) ## 10 samples 
  print(group2) ## 9 samples 
  
  sommut.genes = sort(unique(maf$Hugo_Symbol))
  
  print(sommut.genes) ## 129 genes 
  
  ## --------------------------------------------------------
  
  sommut.stats = NULL
  for(my.gene in sommut.genes) {
    
    print(my.gene)
    my.df = NULL 
    my.stats = NULL 
    
    my.df = unique(maf[maf$Hugo_Symbol %in% my.gene,
                c('Hugo_Symbol','Tumor_Sample_Barcode'), drop=F])
    
    if(nrow(my.df)<2) { next }
    
    rm(x,y)
    x = sum(my.df$Tumor_Sample_Barcode %in% group1)
    y = sum(my.df$Tumor_Sample_Barcode %in% group2)
    
    my.df = matrix(c(x,
                     length(group1) - x, 
                     y,
                     length(group2) - y), nrow = 2)
    
    my.stats = tidy(fisher.test(my.df))
    
    sommut.stats = rbind(sommut.stats,
                         data.frame(Gene = my.gene,
                                    Group1.mut = my.df[1,1],
                                    Group1.nonmut = my.df[2,1],
                                    Group2.mut = my.df[1,2],
                                    Group2.nonmut = my.df[2,2],
                                    my.stats))
    
  }
  dim(sommut.stats)
  
  sommut.stats$p.adj = p.adjust(sommut.stats$p.value, method = 'fdr')
  
  write.csv(sommut.stats,
            file = paste0(maf.file, '.complexVSIsolated.fisher_stats.csv'))
  
}

## --------------------------------------------------------

sessionInfo()
