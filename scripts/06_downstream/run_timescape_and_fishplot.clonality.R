



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

## https://www.bioconductor.org/packages/release/bioc/vignettes/timescape/inst/doc/timescape_vignette.html

# BiocManager::install('timescape')

library(timescape)
library(fishplot)

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
  
  patients.tumor.pdo = c("HN23-10730", "HN24-10810", 
                         "HN24-10856", "HN24-10900", "HN24-10909", "HN24-10936")
  
  print(patients.tumor.pdo)
  
}

## --------------------------------------------

path = 'C:/Users/Owner/OneDrive - University of Pittsburgh/Sunny_work/Projects/TIIL-017-RoseProj/results/exome/PDO_sm61'
setwd(path)

## --------------------------------------------

# ## redirect stdout into local files 
# ## https://stackoverflow.com/questions/7096989/how-to-save-all-console-output-to-file-in-r

# con = NULL 
# con <- file("RETCHER.20250122.log")
# sink(con, append=TRUE)
# sink(con, append=TRUE, type="message")

## --------------------------------------------


my.pt = patients.tumor.pdo[6]
print(my.pt)

# my.pt = 'HN24-10856'

# for(my.pt in patients.tumor.pdo) {
if(TRUE) {
  print(my.pt)
  
  my.indir = NULL 
  
  if(my.pt == 'HN24-10936') {
    
    # my.indir = paste0('pyclone-vi_altDP3/Pt/RETCHER_out_',my.pt,"_1")
    
    my.indir = paste0('pyclone-vi_altDP8/Pt/RETCHER_out/',my.pt,"_7")
    
  } else if (my.pt == 'HN24-10900') {
    
    # my.indir = paste0('pyclone-vi_altDP3/Pt/RETCHER_out/',my.pt)
    # my.indir = paste0('pyclone-vi_altDP3/Pt/RETCHER_out/',my.pt,'_4')
    # my.indir = paste0('pyclone-vi_altDP3/Pt/RETCHER_out/',my.pt,'_16')
    my.indir = paste0('pyclone-vi_altDP3/Pt/RETCHER_out/',my.pt,'_17')
    
    # my.indir = paste0('pyclone-vi_altDP0/Pt/RETCHER_out/',my.pt,'_14')
    # my.indir = paste0('pyclone-vi_altDP0/Pt/RETCHER_out/',my.pt,'_23')
    # my.indir = paste0('pyclone-vi_altDP0/Pt/RETCHER_out/',my.pt,'_15')
    
  } else {
    
    # my.indir = paste0('pyclone-vi_altDP3/Pt/RETCHER_out_',my.pt)
    # my.indir = paste0('pyclone-vi_altDP3/Pt/RETCHER_out/',my.pt,'_16')
    
    # my.indir = paste0('pyclone-vi_altDP0/Pt/RETCHER_out/',my.pt,'_14')
    my.indir = paste0('pyclone-vi_altDP0/Pt/RETCHER_out/',my.pt,'_23')
    
  }
  
  print(my.indir)
  
  ## --------------------------------------------
  
  ## also generate fishplot pdf....
  if(TRUE) {
    
    my.file = NULL 
    my.file = paste0(
      my.indir,'/tree/clone_proportion.tsv'
    )
    
    my.df = NULL 
    my.df = read.delim(my.file)
    
    
    my.df
    
    dim(my.df)
    
    frac.table = NULL 
    frac.table = reshape2::dcast(cluster ~ sample,
                                 value.var = 'mean_ccf',
                                 data = my.df,
                                 fun.aggregate = sum)
    
    row.names(frac.table) = frac.table$cluster
    frac.table = frac.table[,-1]
    
    frac.table = frac.table * 100
    frac.table[frac.table==0] = 0.0001 ## otherwise fishplot errors out if a clone is present in some samples but not others (missing clones); so I just set it as a very small value but >0
    
    ## --------------------------------------------
    
    timepoints = NULL 
    timepoints = 1:ncol(frac.table) ## must be numeric 
    
    ## --------------------------------------------
    
    parents = NULL 
    parents = unique(my.df[,c('cluster','parent')])
    
    ## must be in same order 
    print(all.equal(row.names(frac.table),as.character(parents$cluster)))
    #TRUE
    
    parents = parents$parent ## a vector, not df or matrix
    
    ## --------------------------------------------
    
    fish = NULL 
    
    #create a fish object
    fish = createFishObject(as.matrix(frac.table),parents,timepoints=timepoints)
    
    #calculate the layout of the drawing
    fish = layoutClones(fish)
    
    ## --------------------------------------------
    
    pdf(file = paste0(my.indir, '.fishplot.pdf'),
        height = 5, width = 5)
    
    print(fishPlot(
      fish,
      shape = "spline",
      vlines = timepoints,
      col.vline = "#FFFFFF99",
      vlab = colnames(frac.table),
      border = 0.5,
      col.border = "#777777",
      pad.left = 0.2,
      ramp.angle = 0.5,
      title = NULL,
      title.btm = my.pt,
      cex.title = 0.5,
      cex.vlab = 0.7,
      bg.type = "gradient",
      bg.col = c("bisque", "darkgoldenrod1", "darkorange3")
    ))
    
    dev.off()
  }
  
  ## --------------------------------------------
  
  ## timescape 
  if(TRUE) {
    
    my.file = NULL 
    my.file = paste0(
      my.indir,'/tree/clone_proportion.tsv'
    )
    
    my.df = NULL 
    my.df = read.delim(my.file)
    
    
    my.df
    
    clonal_prev = NULL 
    clonal_prev = data.frame(
      timepoint = as.character(my.df$sample),
      clone_id = as.character(my.df$cluster),
      clonal_prev = as.numeric(as.character(my.df$mean_ccf)) * 100
    )
    
    clonal_prev
    
    ## --------------------------------------------
    
    my.file = NULL 
    my.file = paste0(
      my.indir,'/tree/default_tree.tsv'
    )
    
    my.df = NULL 
    my.df = read.delim(my.file)
    
    
    my.df
    
    tree_edges = NULL 
    tree_edges = data.frame(
      source = as.character(my.df$parent),
      target = as.character(my.df$child)
    )
    
    tree_edges
    
    
    ## --------------------------------------------
    
    my.file = NULL 
    my.file = paste0(
      my.indir,'/cluster/cluster_res.tsv'
    )
    
    my.df = NULL 
    my.df = read.delim(my.file)
    
    
    head(my.df)
    
    ## select TP53 mutations ...
    
    my.df = my.df[grep('[!]TP53[:]',my.df$gene),,drop=F]
    my.df
    # chr     pos                                                 gene
    # 48 chr17 7674221 chr17:7674221:G:A!TP53:p.R248W:Missense_Mutation:SNP
    # TIIL_0165.var TIIL_0167.var TIIL_0168.var TIIL_0190.var TIIL_0260.var
    # 48            85           101           147           107           140
    # TIIL_0165.depth TIIL_0167.depth TIIL_0168.depth TIIL_0190.depth
    # 48             350             102             147             107
    # TIIL_0260.depth TIIL_0165.ccf TIIL_0167.ccf TIIL_0168.ccf TIIL_0190.ccf
    # 48             140        0.7853             1             1             1
    # TIIL_0260.ccf cluster
    # 48             1       1
    
    x = NULL 
    y = NULL 
    x = unique(my.df[,c('chr','pos','gene','cluster'),drop=F])
    y = my.df[,c('gene',colnames(my.df)[grep('^TIIL_',colnames(my.df))])]
    y = reshape2::melt(y, id.vars = c('gene'))
    y$sample_id = gsub('[.]\\S+$','',as.character(y$variable))
    y$metrics = gsub('^\\S+[.]','',as.character(y$variable))
    y = y[y$metrics %in% c('var','depth'),,drop=F]
    y = reshape2::dcast(gene + sample_id ~ metrics,
                        value.var = 'value',
                        data = y,
                        fun.aggregate = sum)
    y$vaf = y$var / y$depth
    
    x$gene = as.character(x$gene)
    y$gene = as.character(y$gene)
    
    x = merge(x, y, by = 'gene')
    
    mutations = data.frame(
      chrom = as.character(x$chr),
      coord = as.numeric(as.character(x$pos)),
      clone_id = as.character(x$cluster),
      timepoint = as.character(x$sample_id),
      VAF = as.numeric(as.character(x$vaf))
      
    )
    
    mutations
    
    ## --------------------------------------------
    
    # pdf(file = paste0('RETCHER_out/',my.pt,'.timescape.pdf'),
    #     height = 5, width = 5)
    
    ## it is interactive!! :) there is a save button on the interactive plot...
    print(timescape(clonal_prev, 
                    tree_edges, 
                    mutations = mutations, ## default = "NA" 
                    clone_colours = "NA",
                    xaxis_title = paste0("Time Point: ", my.pt), 
                    yaxis_title = "Clonal Prevalence",
                    phylogeny_title = "Clonal Phylogeny", 
                    alpha = 50,
                    genotype_position = "centre", ## default = stack
                    perturbations = "NA", 
                    sort = FALSE,
                    show_warnings = TRUE, 
                    width = 900, 
                    height = NULL)
    )
    
    # dev.off()
  }
  
}

# }

## --------------------------------------------

sessionInfo()




