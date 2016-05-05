###################################################################################
##### compute degree and betweenness centrality for a given MITAB2.5 ppi file #####
###################################################################################

get.centrality = function(input.ppi, output.csv){
  # input:
  # - filename of input ppi file, must be a tab-delimited file in MITAB2.5 format
  # - filename of output csv file
  # output:
  # - a .csv file with 4 columns: DIP name of interactor, degree centrality, 
  #   betweenness centrality of the interactor (absolute), and 
  #   betweenness centrality of the interactor (relative) 
  
  ##### read in MITAB2.5 ppi file
  print(paste0('reading in ', input.ppi, '...'))
  ppi = read.delim(input.ppi, header=F, skip=1, as.is=T)
  
  ##### trim 
  ### keep only columns containing info on unique identifiers
  # MITAB2.5 specification avail. at 
  # https://code.google.com/archive/p/psimi/wikis/PsimiTabFormat.wiki
  ppi = ppi[, c(1,2)]
  
  ##### get all unique proteins by DIP name
  # interactor As
  dip.names.a = sapply(ppi[, 1], get.interactor.name)
  names(dip.names.a) = NULL
  ppi[, 1] = dip.names.a
  
  # interactors Bs
  dip.names.b = sapply(ppi[, 2], get.interactor.name)
  names(dip.names.b) = NULL
  ppi[, 2] = dip.names.b
  
  # unique interactors
  dip.uniq = union(dip.names.a, dip.names.b)
  n.uniq = length(dip.uniq)
  
  # num of unique interactors should >= number of unique interactor As or Bs
  stopifnot( (n.uniq >= length(unique(dip.names.a))) & 
             (n.uniq >= length(unique(dip.names.b))) )
  
  ##### build adjacency matrix
  print(paste0('constructing adjacency matrix for ', n.uniq, 
               ' unique interactors (this might take a while)...'))
  adj.mtx = matrix(0, nrow=n.uniq, ncol=n.uniq)
  rownames(adj.mtx) = dip.uniq
  colnames(adj.mtx) = dip.uniq
  
  for (i in 1:nrow(ppi)){
    if (ppi[i, 1]!=ppi[i, 2]){
      # only add edge for non-self interaction
      adj.mtx[ppi[i, 1], ppi[i, 2]] = adj.mtx[ppi[i, 1], ppi[i, 2]] + 1
      adj.mtx[ppi[i, 2], ppi[i, 1]] = adj.mtx[ppi[i, 2], ppi[i, 1]] + 1
    }
  }
  
  # adj.mtx shoule be symmetric and contain no NA
  stopifnot( isSymmetric(adj.mtx) & sum(is.na(adj.mtx))==0 )
  # convert to binary adjacency matrix
  adj.mtx = adj.mtx >= 1
  
  ##### degree centrality
  print('computing degree centrality for each unique interactor...')
  centrality.deg = rowSums(adj.mtx)
  
  ##### betweenness centrality
  print('computing betweenness centrality for each unique interactor...')
  
  ### calculated via Brandes algorithm
  # implemented by brandes.betweenness.centrality() in package 'RBGL
  if ( !('RBGL' %in% installed.packages()[, 'Package']) ){
    # install RBGL from Bioconductor if not installed
    print('package RBGL is required; installing from bioconductor...')
    source("https://bioconductor.org/biocLite.R")
    biocLite("RBGL")
  } else {
    require(RBGL)
  }
  
  ### convert adjacency matrix into a 'graphNEL' object
  adj.graph = as(adj.mtx,"graphNEL")
  
  ### compute betweenness centrality
  brandes = brandes.betweenness.centrality(adj.graph)
  centrality.btw.abs = brandes$betweenness.centrality.vertices
  centrality.btw.rltv = brandes$relative.betweenness.centrality.vertices
  
  ##### output
  print('exporting...')
  centrality = data.frame(interactor = dip.uniq)
  centrality = cbind(centrality, 
                     degree_centrality = as.numeric(centrality.deg),
                     betweenness_centrality = as.numeric(centrality.btw.abs),
                     betweenness_centrality_relative = as.numeric(centrality.btw.rltv))
  
  write.table(centrality, file=output.csv, row.names=F, sep=",", quote=F)
}


###################################################################
##### extract DIP name from unique identifier for interactors #####
###################################################################

get.interactor.name = function(interactor.name){
  # input:
  # - a string containing Unique identifier for interactor A or B
  # output: 
  # - a string containing DIP name extracted from input
  
  # format of interactor name: "DIP-343N|refseq:NP_009971|uniprotkb:P23255"
  # separate by |
  all.names = unlist(strsplit(interactor.name, split="\\|"))
  # get index of DIP name
  dip.idx = base::grepl(pattern='DIP', x=all.names)
  # use only 1 DIP name in case there are multiple (hence the [1])
  dip.name = all.names[dip.idx][1]
  return(dip.name)
}
