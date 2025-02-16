

rm(list = ls())
gc()

## libs
if(TRUE) {
  
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  
  library(maftools)
  
  
}

## ---------------------------------------------

## settings
if(TRUE) {
  
  project = 'PDO_sm61'
  
  ## ---------------------------------------------
  
  patients.tumor.pdo = c("HN23-10730", "HN24-10810", 
                         "HN24-10856", "HN24-10900", "HN24-10909", "HN24-10936")
  
  print(patients.tumor.pdo)
  
  ## ---------------------------------------------
  
  patient.list = readRDS(
    'C:/Users/Owner/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61/PDO_sm61.patient.list.rds'
  )
  
  patient.list
  
}

## ---------------------------------------------

## only run once!
if(FALSE) {
  
  data.path = 'C:/Users/Owner/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61'
  
  ## ---------------------------------------------
  
  ## import altDP8 MAF
  if(TRUE) {
    
    maf.list = list()
    
    for(my.pt in patients.tumor.pdo) {
      
      print(my.pt)
      
      my.file = NULL 
      my.file = paste0(project,'.exome.fixed.gnomAD_exome_AF0.0001_altDP8.maf.',
                       my.pt,'.maf')
      
      my.df = NULL 
      my.df = read.delim(file.path(data.path, my.file))
      
      dim(my.df)
      
      ## ---------------------------------------------
      
      maf.list[[my.pt]] = my.df
    }
    
    names(maf.list)
    
    saveRDS(maf.list,
            file = file.path(data.path,
                             paste0(project,
                                    '.exome.fixed.gnomAD_exome_AF0.0001_altDP8.maf.pt',
                                    length(patients.tumor.pdo),'.maf.list.rds')))
  }
  
  ## ---------------------------------------------
  
  ## import altDP8 pyclone input
  if(TRUE) {
    
    pyclone.in.list = list()
    
    for(my.pt in patients.tumor.pdo) {
      
      print(my.pt)
      
      my.file = NULL 
      my.file = file.path('pyclone-vi_altDP8/Pt/',
                          paste0(project,'.gatkASE.pyclone.list.',
                                 my.pt,'.Pt.txt'))
      
      
      my.df = NULL 
      my.df = read.delim(file.path(data.path, my.file))
      
      dim(my.df)
      
      ## ---------------------------------------------
      
      pyclone.in.list[[my.pt]] = my.df
    }
    
    names(pyclone.in.list)
    
    saveRDS(pyclone.in.list,
            file = file.path(data.path,
                             'pyclone-vi_altDP8/Pt/',
                             paste0(project,
                                    'gatkASE.pyclone.list..pt',
                                    length(patients.tumor.pdo),'.pyclone.in.list.rds')))
    
  }
  
  ## ---------------------------------------------
  
  ## update altDP MAF by rescuing sommuts...
  if(TRUE) {
    

    
  }
  
  
}

## ---------------------------------------------

path = 'C:/Users/Owner/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61'
setwd(path)

## ---------------------------------------------

## data files
if(TRUE) {
  
  pyclone.in.list.file = 'pyclone-vi_altDP8/Pt//PDO_sm61gatkASE.pyclone.list..pt6.pyclone.in.list.rds'
  
  pyclone.in.list = readRDS(pyclone.in.list.file)
  
}

## ---------------------------------------------

## draw plots
if(TRUE) {
  
  data.plot.all = NULL 
  i = 0
  
  for(my.pt in patients.tumor.pdo) {
    i= i+1
    print(my.pt)
    
    my.tumors = NULL 
    my.tumors = patient.list[[my.pt]]
    print(my.tumors)
    
    ## ---------------------------------------------
    
    ## prep inputs
    if(TRUE) {
      
      my.df = NULL 
      my.df = pyclone.in.list[[my.pt]]
      table(my.df$sample_id)
      
      ## split cols 
      head(my.df)
      
      my.df$mutation_id.0 = my.df$mutation_id
      
      my.df = my.df %>%
        tidyr::separate(mutation_id, sep = ':|[!]',
                        c('CHROM','POS','REF','ALT',
                          'Hugo_Symbol','Protein_Change',
                          'Variant_Classification',
                          'Variant_Type'))

      head(my.df)
      
      data.frame(table(my.df$Variant_Classification))
      # Var1 Freq
      # 1             3'UTR  120
      # 2           5'Flank   15
      #       3             5'UTR  135
      # 4            Intron 1620
      # 5 Missense_Mutation 1904
      # 6 Nonsense_Mutation  160
      # 7               RNA   75
      # 8            Silent  845
      # 9       Splice_Site   45
      
      dim(my.df) ## 4919   18
      my.df = my.df[my.df$Variant_Classification %in%
                      c('Missense_Mutation',
                        'Nonsense_Mutation',
                        'Splice_Site'),]
      dim(my.df) ## 2109   18
      
      head(my.df)
      
      ## ---------------------------------------------
      
      ## add flags
      
      my.df$VAF = my.df$alt_counts / 
        (my.df$alt_counts + my.df$ref_counts)
      
      my.df$VAF_presence = 0
      my.df$VAF_presence[which(my.df$VAF > 0)] = 1
      
      table(my.df$Variant_Classification)
      
      my.df$VAF_presence[which(
        my.df$Variant_Classification == 'Missense_Mutation' &
          my.df$VAF_presence >= 1
      )] = 2
      
      my.df$VAF_presence[which(
        my.df$Variant_Classification == 'Nonsense_Mutation'  &
          my.df$VAF_presence >= 1
      )] = -1
      
      my.df$VAF_presence[which(
        my.df$Variant_Classification == 'Splice_Site'  &
          my.df$VAF_presence >= 1
      )] = -2
      
      table(my.df$VAF_presence)
      
      my.df$DP = my.df$ref_counts + my.df$alt_counts
      
      ## ---------------------------------------------
      
      ## which genes to show on the oncoprint?
      
      ## driver genes in HNSCC 
      hnscc.known.genes = read.delim('HNSCC_HPVneg.known_sommut.v2.txt')
      print(hnscc.known.genes)

      
      # ## variants with high VAF in the source tumor 
      # my.genes.select = NULL 
      # my.genes.select = my.df[my.df$sample_id==my.tumors[1],,drop=F]
      # dim(my.genes.select) ## 423  20
      # my.genes.select = my.genes.select[my.genes.select$DP >= 10 &
      #                           my.genes.select$alt_counts >= 3,,drop=F]
      # dim(my.genes.select) ## 234  20
      # 
      # my.genes.select = my.genes.select[rev(order(my.genes.select$VAF)),,drop=F]
      # 
      # my.genes.select = my.genes.select[!my.genes.select$Hugo_Symbol %in%
      #                                     hnscc.known.genes$Gene,,drop=F]
      # 
      # head(my.genes.select)
      
      my.df = my.df[my.df$Hugo_Symbol %in% hnscc.known.genes$Gene,,drop=F]
      dim(my.df) ## 75 20
      
      ## ---------------------------------------------
      
      
      data.plot = NULL 
      data.plot = reshape2::dcast(
        mutation_id.0 ~ sample_id,
        value.var = 'VAF_presence',
        data = my.df,
        fun.aggregate = sum
      )
      
      head(data.plot)

      row.names(data.plot) = data.plot$mutation_id.0
      data.plot = data.plot[,-1]
      
      row.names(data.plot) = 
        gsub('^\\S+[!]','',row.names(data.plot))
      
      row.names(data.plot) = gsub(':SNP','',row.names(data.plot))
      
      row.names(data.plot) = gsub('[:][A-Z][A-Za-z_]+','',row.names(data.plot))
      
      row.names(data.plot) = gsub('[:][.]$',':Spice_Site',row.names(data.plot))
      
      head(row.names(data.plot))
      
      if(i == 1) {
        data.plot.all = data.plot
      } else {
        
        data.plot.all = merge(
          data.plot.all, data.plot, 
          by = 'row.names',
          all = T
        )
        
        row.names(data.plot.all) = data.plot.all[,1]
        data.plot.all = data.plot.all[,-1]
      }

    }
    
    ## ---------------------------------------------
    
    ## heatmap
    if(TRUE) {
      
      library(circlize)
      library(ComplexHeatmap)
      library(plotrix)
      library(RColorBrewer)
      
      centered_data = NULL 
      centered_data = data.plot
      
      myheatcol = colorRamp2(c(-2, 0, 2), c("blue","white", "red"))
      col.title = paste0(my.pt)
      row.title = paste0('Somatic Mutations')
      
      p1 = NULL 
      p1 = Heatmap(centered_data,
                   na_col = "#000000",
                   col = myheatcol,
                   rect_gp = gpar(col = '#ffffff'),
                   show_heatmap_legend = T,
                   column_title = col.title,
                   row_title = row.title,
                   # column_title_side = 'bottom',
                   column_names_side = 'top',
                   row_dend_width = unit(5, "cm"),
                   column_dend_height = unit(5, "cm"),
                   # km = 2,
                   cluster_rows = F,
                   cluster_columns = F,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "ward.D2",
                   clustering_distance_columns = "euclidean",
                   clustering_method_columns = "ward.D2",
                   show_row_names = T,
                   show_column_names = T,
                   # top_annotation = plot.anno,
                   heatmap_legend_param = list(title = 'Mutations', 
                                               color_bar = "discrete")
      )
      
      p1
      
      pdf(file = paste0(my.pt,'.sommut.heatmap.pdf'), width = 5, height = 5)
      print(p1)
      dev.off()
      
      
    }
  }
  
  ## ---------------------------------------------
  
  head(data.plot.all)
  data.plot.all[is.na(data.plot.all)] = 0
  
  plot.order = NULL 
  plot.order = data.frame(
    mutation_id = row.names(data.plot.all),
    stringsAsFactors = F)
  
  plot.order$Hugo_Symbol = gsub('[:]\\S+$','',plot.order$mutation_id)
  
  x = NULL 
  for(my.gene in hnscc.known.genes$Gene) {
    
    x = c(x,
          plot.order$mutation_id[plot.order$Hugo_Symbol==my.gene])
  }
  
  x
  
  plot.order = plot.order[order(match(plot.order$mutation_id,x)),,drop=F]
  
  head(plot.order)
  
  data.plot.all = data.plot.all[order(match(
    row.names(data.plot.all), plot.order$mutation_id
  )),,drop=F]
  
  ## ---------------------------------------------
  
  ## heatmap
  if(TRUE) {
    
    library(circlize)
    library(ComplexHeatmap)
    library(plotrix)
    library(RColorBrewer)
    
    centered_data = NULL 
    centered_data = data.plot.all
    
    myheatcol = colorRamp2(c(-2, 0, 2), c("blue","white", "red"))
    col.title = paste0(my.project)
    row.title = paste0('Somatic Mutations')
    
    p1 = NULL 
    p1 = Heatmap(centered_data,
                 na_col = "#000000",
                 col = myheatcol,
                 rect_gp = gpar(col = '#ffffff'),
                 show_heatmap_legend = T,
                 column_title = col.title,
                 row_title = row.title,
                 # column_title_side = 'bottom',
                 column_names_side = 'top',
                 row_dend_width = unit(5, "cm"),
                 column_dend_height = unit(5, "cm"),
                 # km = 2,
                 cluster_rows = F,
                 cluster_columns = F,
                 clustering_distance_rows = "euclidean",
                 clustering_method_rows = "ward.D2",
                 clustering_distance_columns = "euclidean",
                 clustering_method_columns = "ward.D2",
                 show_row_names = T,
                 show_column_names = T,
                 # top_annotation = plot.anno,
                 heatmap_legend_param = list(title = 'Mutations', 
                                             color_bar = "discrete")
    )
    
    p1
    
    pdf(file = paste0(my.project,'.sommut.heatmap.pdf'), width = 9, height = 6)
    print(p1)
    dev.off()
    
    
  }
  
}

## ---------------------------------------------


