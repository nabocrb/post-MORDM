### SOM grid search experiment
### Nathan Bonham
### 8/30/21

library(kohonen) # SOM
library(aweSOM) # quality metrics: quantization error (perc. variance), topo error
library(lhs) # latin hypercube design of SOM parameters
library(clusterSim) # Davies-Bouldin index
rm(list=ls())

source("G:/My Drive/CU Boulder/CRB publications/Uncertainty characterization and RDM sensitivity paper/R/Library.R")

# Import DV
setwd('G:/My Drive/CU Boulder/CRB publications/Uncertainty characterization and RDM sensitivity paper/R/data')
#setwd('G:/My Drive/CU Boulder/CRB publications/Uncertainty characterization and RDM sensitivity paper/R/data')
Archive='Archive_463_Condensed.txt'
DV=read.table(Archive, header = T, sep = "")

#### k means clustering to get a sense for how many nodes we might want
# ALL POLICIES

scale.optimization=scale(optimization[-1])

DBindex=2:100

set.seed(7)

for(k in 2:100){
  
  temp=kmeans(x = scale.optimization, centers = k, nstart = 25)
  cluster=temp$cluster
  temp=index.DB(x = scale.optimization,cl = cluster)  
  DBindex[k-1]=temp$DB
}

DB.df=data.frame(k=2:100, DBindex)

ggplot(data = DB.df, aes(x=k, y=DBindex))+
  geom_line()+
  ggtitle("Davies-Bouldin index for optimization objectives")+
  ylab("index")+
  xlab("clusters")+
  scale_x_continuous(breaks=seq(0,100, by=10))+
  geom_point(x=13, y=DB.df$DBindex[12], size=3, color="blue")+
  geom_text(x=15, y=DB.df$DBindex[12]-.01, label="k=13")

# conclusion: Davies Bouldin is lowest at k = 13 except where k > 80.
# Since the index values are not much better above k > 80, let's focus closer to k=13.
# Let's search from 4 to 25 nodes with SOM and look at k=13 afterwards

######################################################################
########################## lhs #######################################

n_samples=1000

# testing radius, grid size, distance function, neighborhood function, toroidal or planar

set.seed(7)

#unit_cube=geneticLHS(n=n_samples, k=4, criterium = 'Maximin') # optimizes Maximin optimality criterion.s
unit_cube=improvedLHS(n=n_samples, k=5) # optimizes for euclidean distance between points

unit_cube=data.frame(unit_cube)
colnames(unit_cube)=c("radius", "n_nodes", "distance_fnc", "neighbor_fnc", "toroidal")
# convert to uniform sample across given ranges

radius=c(0,1) # as a fraction of max unit-to-unit distances
n_nodes=c(4, 25)
distance_fnc=c("sumofsquares", "euclidean", "manhattan")
neighbor_fnc=c("bubble", "gaussian")
toroidal=c(T, F)
cube=unit_cube

cube$radius=runif(unit_cube$radius, min=radius[1], max=radius[2])
cube$n_nodes=round(runif(unit_cube$n_nodes, min=n_nodes[1], max=n_nodes[2]))

### transform distance functions to discrete
probs=seq(0,1, length.out = (length(distance_fnc)+1)) # sets equal probability ranges for finding which 0 to 1 values correspond to given distance function

for (i in 1:(length(probs)-1)){
  
  cube$distance_fnc[which(unit_cube$distance_fnc > probs[i] & unit_cube$distance_fnc<= probs[i+1])]=distance_fnc[i]
    
    
}

### transform neighbor functions and toroidal TF to discrete

probs=seq(0,1, length.out = (length(neighbor_fnc)+1)) # sets equal probability ranges for finding which 0 to 1 values correspond to given neighbor function

for (i in 1:(length(probs)-1)){
  
  cube$neighbor_fnc[which(unit_cube$neighbor_fnc > probs[i] & unit_cube$neighbor_fnc<= probs[i+1])]=neighbor_fnc[i]
  cube$toroidal[which(unit_cube$toroidal > probs[i] & unit_cube$toroidal<= probs[i+1])]=toroidal[i]
  
  
}

## test that each value is equally represented in the hypercube

fraction=1:length(distance_fnc) # preallocate
c=1
for (i in distance_fnc){
  fraction[c]=sum(cube$distance_fnc==i)/nrow(cube)
  c=c+1
}

fraction # each element in fraction should be equal. It is the fraction of your samples that correspond to each demand value

fraction=1:length(neighbor_fnc) # preallocate
c=1
for (i in neighbor_fnc){
  fraction[c]=sum(cube$neighbor_fnc==i)/nrow(cube)
  c=c+1
}

fraction

fraction=1:length(toroidal) # preallocate
c=1
for (i in toroidal){
  fraction[c]=sum(cube$toroidal==i)/nrow(cube)
  c=c+1
}

fraction

##################################################################################################################
############################## calculate SOM map and report quality metrics for each #############################

scaled_data= scale.optimization

# Kohonen 2013 - The Essentials of the SOm - section 3.5. 
# "It is advisable to select the lengths of the horizontal and vertical dimensions of the array
# to correspond to the lengths of the two largest principal components (ie those with the highest 
# eigenvalues of the input correlation matrix)". So, I calculate the covariance matrix, which is square,
# then compute the eigenvalues, then take ratio of top two
temp.cov=cov(scaled_data)
temp.eigenValues=eigen(temp.cov)$values
topEV=sort(temp.eigenValues,decreasing = T)[1:2]

ratio=topEV[1]/topEV[2]

# create initialization matrix as uniformly sampled plane along first two PCs. See Kohonen 2013 and Clark et al 2013

X=scaled_data
temp.eigenVectors=eigen(cov(X))
RM=temp.eigenVectors$vectors[,1:2] # rotation matrix, which is the eigenvectors corresponding to first and second largest eigenvalues
PCs=X %*% RM
PC1range=range(PCs[,1]) # will use these to sample uniformly along PC1 and PC2 within the for loop
PC2range=range(PCs[,2])

metric.df=matrix(nrow=n_samples, ncol= 4) # preallocate matrix to store quality of fit metrics
metric.df=data.frame(metric.df)

params.df=matrix(nrow=n_samples, ncol=6)
params.df=data.frame(params.df)
colnames(params.df)=c("x_dim", "y_dim", "radius", "dist_fnc", "neighbor_fnc", "toroidal")




for (i in 1:nrow(cube)){
  
  # calculate x and y using the ratio calculated above: 
  # eqn 1: N =x*y, where N is n_nodes
  # eqn2: x/y = e1/e2 = ratio
  # solving for x:
  
  N=cube$n_nodes[i]
  
  x=sqrt(N*ratio)
  y=N/x
  
  x=round(x)
  y=round(y)

  
  # create nueron initialization matrix
  
  D1=seq(PC1range[1], PC1range[2], length.out = x) # sequence of nuerons along PC1
  D2=seq(PC2range[1], PC2range[2], length.out = y) # sequence of nuerons along PC1
  
  PG=expand.grid(D1,D2) # create rectangular matrix where D1 is repeated for every value of D2, projected in PC space
  
  # plot(x=PG$Var1,y=PG$Var2, type = "p" ) # to see example of the grid in PC space
  # points(PCs[,1], PCs[,2], col="red") # to see the policies in PC space
  # unproject the projected grid (PG) back into original data space
  
  IG=as.matrix(PG) %*% t(RM) # neuron initialization matrix
  
  # fit SOM
  
  radius=fraction2radius(fraction = cube$radius[i], x = x, y = y, shape = "hexagonal")

  temp.SOM=som(X =scaled_data,radius= radius, dist.fcts=cube$distance_fnc[i],
               grid=somgrid(xdim=x, ydim=y, topo = "hexagonal", toroidal = cube$toroidal[i], neighbourhood.fct = cube$neighbor_fnc[i]),
               rlen=10000, keep.data=T, init=IG, mode="batch") # Kohonen 2013 recommends the batch algorithm vs the original stepwise recursive approach "because it is faster and safer".

  # test number of training iterations needed
  
  #plot(temp.SOM, type = "changes")
  
  temp.metrics=somQuality(temp.SOM, traindat = scaled_data)
  
  #temp.metrics[1:4] to print to console if testing something
  
  metric.df[i,]=temp.metrics[1:4]
  
  if(i==1){
    colnames(metric.df)=names(temp.metrics)[1:4]
  }
  
  params.df[i,]=c(x,y,radius, cube$distance_fnc[i], cube$neighbor_fnc[i], cube$toroidal[i])
  
}

setwd('G:/My Drive/CU Boulder/CRB publications/Uncertainty characterization and RDM sensitivity paper/R/results/Optimization - all policies')
write.table(metric.df, file="SOM quality metrics.txt")
write.table(params.df, file = "SOM parameters.txt")
write.table(cube, file="Hypercube design.txt")






















