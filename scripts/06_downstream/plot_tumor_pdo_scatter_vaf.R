
rm(list = ls())
gc()

## libs 
if(TRUE) {
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  
  library(dplyr)
}

## -------------------------------------------

## settings
if(TRUE) {
  patients.tumor.pdo = c(
    'HN23-10730','HN24-10810','HN24-10856','HN24-10900','HN24-10909','HN24-10936'
  )
  
  patients.tumor.pdo.cluster.runs = list(
    `HN23-10730` = c('altDP0','HN23-10730_23'),
    `HN24-10810` = c('altDP0','HN24-10810_23'),
    `HN24-10856` = c('altDP0','HN24-10856_23'),
    `HN24-10900` = c('altDP3','HN24-10900_17'),
    `HN24-10909` = c('altDP0','HN24-10909_23'),
    `HN24-10936` = c('altDP8','HN24-10936_7')
  )

}

## -------------------------------------------

path = 'C:/Users/Owner/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61'
setwd(path)

## -------------------------------------------

# my.pt = patients.tumor.pdo[6]
# print(my.pt)

for(my.pt in patients.tumor.pdo) {
# if(TRUE) {
  
  print(my.pt)
  
  ## -------------------------------------------
  
  ## import data
  if(TRUE) {
    
    my.file = NULL 
    
    ## #
    my.file = paste0(
      'pyclone-vi_',
      ## use most st# or, use the same variant files for RETCHER runringent cutoff for variants...
      # 'altDP8', 
      ## or, use the same variant files for RETCHER run
      patients.tumor.pdo.cluster.runs[[my.pt]][1], 
      '/Pt/',
      'PDO_sm61.gatkASE.pyclone.list.',my.pt,'.Pt.txt'
    )
    print(my.file)

    my.df = NULL 
    my.df = read.delim(my.file)
    dim(my.df) ## 4919    9
    
    table(my.df$sample_id)
    # IIL_0285 TIIL_0287 TIIL_0288 TIIL_0289 TIIL_0290 
    # 985       983       983       985       983 
    
    my.df$alt.frac = my.df$alt_counts / (
      my.df$alt_counts + my.df$ref_counts
    )
    
    ## -------------------------------------------
    
    my.df.list = list()
    my.tumors = sort(unique(my.df$sample_id))
    print(my.tumors)
    
    ## label low 
    for(my.sm in my.tumors) {
      
      my.df.2 = NULL 
      my.df.2 = my.df[my.df$sample_id==my.sm,,drop=F]
      my.df.2$alt_counts_low = 'high'
      my.df.2$alt_counts_low[which(my.df.2$alt_counts < 8)] = 'low'
      
      my.df.list[[my.sm]] = my.df.2
      
    }
    
    print(names(my.df.list))

    ## -------------------------------------------
    
    ## add clusters!
    ## RETCHER clusters: I spent a lot of picking out parameters ......
    my.cluster.file = NULL 
    my.cluster.file = paste0(
      'pyclone-vi_',
      patients.tumor.pdo.cluster.runs[[my.pt]][1],'/',
      'Pt/RETCHER_out/',
      patients.tumor.pdo.cluster.runs[[my.pt]][2],
      '/cluster/cluster_res.tsv'
    )
    my.cluster = NULL 
    my.cluster = read.delim(my.cluster.file)
    dim(my.cluster) ## 741  19
    table(my.cluster$cluster)
    # 1   2 
    # 76 665
    
    my.cluster$mutation_id = as.character(my.cluster$gene)
    
  }
  
  ## -------------------------------------------
  
  ## prep inputs
  if(TRUE) {
    
    data.plot = NULL 
    data.plot = reshape2::dcast(mutation_id ~ sample_id,
                                value.var = 'alt.frac',
                                data = my.df,
                                fun.aggregate = sum)
    head(data.plot)
    dim(data.plot) ## 985   6

    ## -------------------------------------------
    
    data.plot$mutation_id = as.character(data.plot$mutation_id)
    
    dim(data.plot) ## 985   6
    data.plot = merge(data.plot,
                      my.cluster[,c('mutation_id','cluster')],
                      by = 'mutation_id',
                      all.x = T)
    dim(data.plot) ## 985   7
    
    head(data.plot)
    data.plot$cluster[is.na(data.plot$cluster)] = 99 ## excluded during RETCHER analysis
    
    data.plot$cluster = paste0('Clr',data.plot$cluster)
    
    ## add low count flag in source tumor ...
    data.plot = merge(data.plot,
                      my.df.list[[my.tumors[1]]][,c('mutation_id',
                                                   'ref_counts',
                                                   'alt_counts',
                                                   'alt_counts_low')],
                      by = 'mutation_id',
                      all.x = T)
    dim(data.plot) ## 985   8
    sum(is.na(data.plot$alt_counts_low)) ## 0; all in!
    
    write.csv(data.plot,
              file = paste0(my.file, '.data.csv'))
    
  }
    
  ## -------------------------------------------
  
  ## draw plots
  if(TRUE) {
    
    print(my.tumors)

    data.plot.2 = NULL 
    data.plot.2 = data.plot[!data.plot$cluster %in% c('Clr99'),]
    
    ## -------------------------------------------
    
    p0 = ggplot(data.plot.2, aes(.data[[my.tumors[1]]], 
                                 .data[[rev(my.tumors)[2]]])) +
      geom_point(aes(fill = alt_counts_low,
                     size = ref_counts),
                 shape = 21, alpha = 0.5) + 
      ggtitle(my.pt) +
      xlim(0,1)+
      ylim(0,1) +
      theme_pubr()
    
    p0 ## not necessarily low total DP are the sites with 0 frac in source tumor??
    
    pdf(file = paste0(my.file, '.scatter.01.pdf'),
        width = 5, height = 5)
    print(p0)
    dev.off()
    
    ## -------------------------------------------
    
    p1 = ggplot(data.plot.2, aes(.data[[my.tumors[1]]], 
                                 .data[[rev(my.tumors)[2]]])) +
      geom_point(aes(fill = cluster),
                 shape = 21, alpha = 0.5) + 
      ggtitle(my.pt) +
      stat_cor(method = 'spearman') +
      xlim(0,1)+
      ylim(0,1)+
      theme_pubr()
    
    
    p2 = ggplot(data.plot.2, aes(.data[[my.tumors[1]]], 
                                 .data[[rev(my.tumors)[1]]])) +
      geom_point(aes(fill = cluster,
                     shape = cluster),
                 shape = 21, alpha = 0.5) +  
      ggtitle(my.pt)  +
      stat_cor(method = 'spearman') +
      xlim(0,1)+
      ylim(0,1)+
      theme_pubr()
    
    p3 = ggplot(data.plot.2, aes(.data[[rev(my.tumors)[2]]],
                                 .data[[rev(my.tumors)[1]]])) +
      geom_point(aes(fill = cluster),
                 shape = 21, alpha = 0.5) + 
      ggtitle(my.pt) +
      stat_cor(method = 'spearman') +
      xlim(0,1)+
      ylim(0,1)+
      theme_pubr()
    
    p1 + p2 + p3
    
    pdf(file = paste0(my.file, '.scatter.02.pdf'),
        width = 15, height = 5)
    print(p1 + p2 + p3)
    dev.off()
  }
  
  

}



