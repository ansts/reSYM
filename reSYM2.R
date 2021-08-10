require(stringdist)
require(Biostrings)
require(parallel)
require(matrixStats)
require(igraph)
require(reshape2)



##########################################################
# ODX  graph
##########################################################
load(file="ODxs_outord_r") # complete results from Gibb's clustering
ODXAdjL=lapply(ODxs_outord,function(m){
  x=m$BGscore
  n=m$BGgroup+1
  y=aggregate(x, by=list(n), FUN=length)
  z=y[y[,2]>8,2]
  names(z)=y[y[,2]>8,1]
  return(z)
})
names(ODXAdjL)=1:790
ODXAdjLw=lapply(ODXAdjL,function(l){
  names(l)=NULL
  return(l)
})
ODXAdjL=lapply(ODXAdjL,function(l){
  return(names(l))
})
ODXel=melt(ODXAdjL)
ODXel=as.matrix(ODXel)
ODXw=(unlist(melt(ODXAdjLw)[,1]))
GODX=graph_from_edgelist(ODXel, directed = F)
E(GODX)$weight=ODXw # score is used as edge weight
GODX=simplify(GODX)
write.graph(GODX,file="GODX.graphml", format="graphml")
GODX=delete.edges(GODX, which(E(GODX)$weight<34)) 

# do it in Gephi
# gephires=read.csv(file="GODX.nodesxx.csv") # number of connected components = 727 at edge weight >33.846
# geres=gephires[,4:6]
# the727L=aggregate(geres[,1],by=list(geres[,2]), FUN=list)

# or in R

geres=components(GODX)$membership
the727L=aggregate(names(geres),by=list(geres), FUN=list)
the727L=the727L[order(lengths(the727L[,2]), decreasing = T),]
the727=the727L[,2]
the727[[727]]=the727[[1]][c(1,2,5)]
the727[[1]]=the727[[1]][-c(1,2,5)]
j=as.numeric(unlist(the727[lengths(the727)==1]))-1
newchosen=sapply(ODxs_outord[j],function(l){
  l[1,2]
})
j=which(lengths(the727)>1)

x=sapply(j,function(j){
  x=the727[[j]]
  p=ODXs_790_30[ODXs_790_30[,2] %in% x,1]
  pAdjL=adjL(p)
  pel=melt(pAdjL)
  pel=as.matrix(pel)
  Gp=graph_from_edgelist(pel, directed = F)
  Gp=simplify(Gp)
  eCGp=eigen_centrality(Gp)
  x=names(eCGp$vector)[which.max(eCGp$vector)]
  return(x)
})
newchosen=c(x,newchosen)
ncAdjL=adjL(newchosen) # check similarities of selected sequences
ncel=melt(ncAdjL)
ncel=as.matrix(ncel)
Gnc=graph_from_edgelist(ncel, directed = F)
Gnc=simplify(Gnc)
compGnc=components(Gnc)
names(compGnc$membership[compGnc$membership>2])
write.csv(newchosen, file="newchosen.csv")
# HCCWTKM and QLSMSSR are excluded to reach 724