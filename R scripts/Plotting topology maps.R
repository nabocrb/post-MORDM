### post-MORDM: figures in main text and appendix
### coded by Nathan Bonham
### March 2022

# load packages and custom functions
source(here("R scripts","Library.R")) 

###############################################################################################
######################## Plot SOM fit metrics, choose parameter set ####################################

cube=read.table(here("case study data","Hypercube design.txt")) # the latin hypercube design of SOM parameters
SOMparameters=read.table(here("case study data","SOM parameters.txt")) # SOM parameters as exactly inputted to SOM. Some parameters are calculated and/or rounded from cube
SOMparameters$nodes=SOMparameters$x_dim*SOMparameters$y_dim
SOMfit=read.table(here("case study data","SOM quality metrics.txt")) # goodness of SOM fit, calculated for every sample using aweSOM package

SOMfit_plot=data.frame(SOMparameters, SOMfit)
colnames(SOMfit_plot)=c("x", "y", "radius", "distance.fnc","neighbor.fnc", "toroidal", "nodes", "quant.e", "perc.variance", "topo.e", "KL.e")
SOMfit_plot$radius=cube$radius

# assign integer to discrete values for plotting
SOMfit_plot$distance.fnc[SOMfit_plot$distance.fnc=="euclidean"]=1
SOMfit_plot$distance.fnc[SOMfit_plot$distance.fnc=="manhattan"]=2
SOMfit_plot$distance.fnc[SOMfit_plot$distance.fnc=="sumofsquares"]=3

SOMfit_plot$neighbor.fnc[SOMfit_plot$neighbor.fnc=="bubble"]=1
SOMfit_plot$neighbor.fnc[SOMfit_plot$neighbor.fnc=="gaussian"]=2

SOMfit_plot$distance.fnc=as.numeric(SOMfit_plot$distance.fnc)
SOMfit_plot$neighbor.fnc=as.numeric(SOMfit_plot$neighbor.fnc)

str(SOMfit_plot)

cat_names=list()
cat_names[[1]]=c("euclidean", "manhattan", "sumofsquares")
cat_names[[2]]=c("bubble", "gaussian")

temp.df=data.frame(SOMfit_plot$nodes, SOMfit_plot$quant.e, SOMfit_plot$topo.e)

# remove dominated parameter sets
dominated=ecr::dominated(t(as.matrix(temp.df))) # returns TRUE if dominated, false if non-dominated

SOMfitFront1=SOMfit_plot[dominated==F,]

# Appendix A5
par_coords(data = SOMfitFront1[,-c(8,11)], max_cols = which(colnames(SOMfitFront1[,-c(8,11)])=="perc.variance"), n_var = 9,
           color_var = "nodes", title = "SOM: parameters and fit metrics", categorical_cols = which(colnames(SOMfitFront1) %in% c("distance.fnc", "neighbor.fnc")  ),
           categorical_names_list = cat_names, legendTF = T, colorbarTitle = "Nodes", img_width = 800, color_scale = "Portland")

### select SOM from shopping with Par Coords

x=5; y=3; distance.fnc="sumofsquares"; toroidal=F; neighbor.fnc="bubble"

radius=SOMfitFront1$radius[which(SOMfitFront1$nodes==15 & SOMfitFront1$quant.e <0.8 &
                                   SOMfitFront1$topo.e<0.03 & SOMfitFront1$distance.fnc==3 &
                                   SOMfitFront1$neighbor.fnc==1)][1]
save_radius=radius

# 91.02 percent variance explained, 0.0259 topo.e

my_vals=SOMfitFront1[which(SOMfitFront1$nodes==15 & SOMfitFront1$quant.e <0.8 &
                             SOMfitFront1$topo.e<0.03 & SOMfitFront1$distance.fnc==3 &
                             SOMfitFront1$neighbor.fnc==1)[1],]

n_misplaced=my_vals$topo.e*nrow(optimization) # 12 policies constitute a topo error

##################################################################################################
################################### create SOM object based on chosen parameter set ##############

scaled_data= as.matrix(scale(optimization[,-1]))

# create initialization matrix as uniformly sampled plane along first two PCs. See Kohonen 2013 and Clark et al 2013

X=scaled_data
temp.eigenVectors=eigen(cov(X))
RM=temp.eigenVectors$vectors[,1:2] # rotation matrix, which is the eigenvectors corresponding to first and second largest eigenvalues
PCs=X %*% RM
PC1range=range(PCs[,1]) # will use these to sample uniformly along PC1 and PC2 within the for loop
PC2range=range(PCs[,2])
D1=seq(PC1range[1], PC1range[2], length.out = x) # sequence of nuerons along PC1
D2=seq(PC2range[1], PC2range[2], length.out = y) # sequence of nuerons along PC1

PG=expand.grid(D1,D2) # create rectangular matrix where D1 is repeated for every value of D2, projected in PC space

plot(x=PG$Var1,y=PG$Var2, type = "p" ) # to see example of the grid in PC space
points(PCs[,1], PCs[,2], col="red") # to see the policies in PC space
# unproject the projected grid (PG) back into original data space

IG=as.matrix(PG) %*% t(RM) # neuron initialization matrix

# sort and rank eigen vectors to understand the x and y axis of SOM
RM=data.frame(RM)
RM$objective=c("M1000", "LB.Dur", "LB.Freq", "LB.Avg", "LB.Max", "P3490", "P.WYR", "LF.Deficit")
E1=dplyr::arrange(RM, -abs(X1)) # sort largest to smallest absolute eigenvalue
E1=E1[,-2] # remove 2nd eigen vector
E1$loading=E1$X1*sqrt(temp.eigenVectors$values[1]) # loading score = eigenvector * sqrt eigenvalue

E2=dplyr::arrange(RM, -abs(X2)) # sort largest to smallest absolute eigenvalue
E2=E2[,-1] # remove 1st eigen vector
E2$loading=E2$X2*sqrt(temp.eigenVectors$values[2]) # loading score = eigenvector * sqrt eigenvalue

# fit SOM

radius=fraction2radius(fraction = radius , x = x, y = y, shape = "hexagonal")

SOM.opt=som(X =scaled_data,radius= radius, dist.fcts=distance.fnc,
    grid=somgrid(xdim=x, ydim=y, topo = "hexagonal", toroidal = toroidal, neighbourhood.fct = neighbor.fnc), rlen=10000, keep.data=T,
    mode="batch", init=IG)

temp=aweSOM::somQuality(SOM.opt, traindat = scaled_data)
firstSOMfit=c("optimization",temp$err.quant, temp$err.varratio, temp$err.topo, "SOM fit")

################################################################################
##################### Plotting topology maps #######################################

### radar plots

scaled_opt=apply(optimization[,-1], MARGIN = 2, FUN = scale1) # scale data 0-1
scaled_opt=data.frame(scaled_opt)
scaled_opt$node=SOM.opt$unit.classif

colnames(scaled_opt)=c("M1000","LB.Dur","LB.Freq","LB.Avg","LB.Max","P3490","P.WYR", "LF.Deficit", "node")

my_order=c("P.WYR", "LF.Deficit", "LB.Dur", "LB.Max", "LB.Freq","LB.Avg", "M1000", "P3490", "node")

scaled_opt=scaled_opt[,my_order]


## plot Figure 5

radar_grid(x=x, y=y, data=scaled_opt, color = "lightblue", line_color = "black", alpha=0.3, hex_shift = "right") # programatically create hexagon grid
# radar_grid(x=x, y=y, data=scaled_opt, color = "lightblue", line_color = "black", alpha=0.3, hex_shift = F) # place neurons in rectangular grid, edit with inkscape

mtext(text="objectives", outer=T, cex=1.3, line = -1.5,adj = 0.415)

#### component planes, Appendix A6

filter.optimization=optimization[,-1] # also remove ID column
filter.optimization$node=SOM.opt$unit.classif
colnames(filter.optimization)=c("Mead.1000","LB.Dur","LB.Freq","LB.Avg","LB.Max","Powell.3490","Powell.WYR", "LF.Deficit", "node")
my_order=c("LF.Deficit","Powell.WYR", "Powell.3490", "Mead.1000","LB.Dur","LB.Avg", "LB.Max", "LB.Freq")
filter.optimization=filter.optimization[,c(my_order, "node")]

units=c("%", "MAF", "%", "%", "Years","KAF", "KAF", "%")
optimization_circles=SOM_bubbles(SOMlist = SOM.opt ,data=filter.optimization, units = units, my_order = my_order, title = "Optimization objectives",ncol = 2, nrow = 4, hex=T, labels=T)
optimization_circles

##########################################################################################
############### robustness calculations and plotting relationships #######################

# note: robustness calculations are used in main text and appendix

# clean up robustness simulations data with shorter names

colnames(obj_all)[-c(1,2)]=c("LB.Avg","LB.Avg.Policy","LF.Deficit","M1000","P3490","LB.Dur","LB.Freq", "LB.Max", "P.WYR")
obj_all=obj_all[-which(colnames(obj_all)=="LB.Avg.Policy")]
my_order=c("TraceNumber","policy","LF.Deficit","P.WYR", "P3490", "M1000","LB.Dur","LB.Avg", "LB.Max", "LB.Freq")
obj_all=obj_all[,my_order]

obj_all$P.WYR=obj_all$P.WYR/1e6 # convert AF to MAF
obj_all[,c("LB.Avg", "LB.Max")]=obj_all[,c("LB.Avg", "LB.Max")]/1e3 # convert AF to KAF

#### Testing multiple robustness metrics

robustness.list=list()
my_objs=colnames(obj_all)[-c(1,2)]

### Laplace's Principle of Insufficient Reason (mean)

robustness.list[["mean"]]=LaplacePIR(data = obj_all, objectives = my_objs)

### Maximin (worst case scenario, try 90th and 100th percentile)

robustness.list[["90% maximin"]]=maximin(data= obj_all, objectives = my_objs, percentile = rep(90, ncol(obj_all)-2) )
robustness.list[["maximin"]]=maximin(data= obj_all, objectives = my_objs, percentile = rep(100, ncol(obj_all)-2) )

### Regret from best, no normalization of regret 

regret.df=data.frame(matrix(nrow=nrow(optimization), ncol = length(my_objs)+1))
colnames(regret.df)=c("policy",my_objs)
regret.df$policy=1:nrow(optimization)
for (i in my_objs){
  
  regret.df[,i]=regret2(data = obj_all, objectives = i, SOW_agg_method = "percentile",
                        percentile = 0.9, best_if = rep("min", length(my_objs)), Obj_agg_method = "none" , scale = F)$regret2[,2]
  
}

robustness.list[["90% regret from best"]]=regret.df


#### decision maker scenario

## Water delivery

# satisficing if LB.Avg<=600 KAF and LB.Dur <= 10 years

sat.LBAvg.LBDur=satisficing(data=obj_all, objectives = c("LB.Avg", "LB.Dur"), thresholds = c(600, 10), fail_if_inequality = c("greater", "greater"))

## Reservoir storage

# satisficing if both M1000 <= 10% and P3490 <= 5%

sat.P3490.M1000=satisficing(data = obj_all, objectives = c("P3490", "M1000"), thresholds = c(5,10), fail_if_inequality = c("greater", "greater"))

robustness.list[["Stakeholder satisficing metrics"]]=data.frame(Delivery=sat.LBAvg.LBDur$satisficing, Storage=sat.P3490.M1000$satisficing)

##############################################################################################
########################### topology maps of robustness metrics ##############################

### component planes (Figure 7, Appendix A11-13 part b)

robustness.planes=list()

stakeholder_units=c("fraction SOW", "fraction SOW")
stakeholder_min=c(F, F)

for(i in names(robustness.list)){
  
  w=14.75; h=5.3 # set width and height of figures
  
  robustness.df=robustness.list[[i]]
  
  robustness.df$node=SOM.opt$unit.classif
  
  if(ncol(robustness.df)==10){
    
    robustness.df=robustness.df[,-1] # remove policy column
    
    my_units=c("%", "MAF", "%", "%", "years", "KAF", "KAF", "%")
    cols=4
    rows=2
    min_cols=rep(T, length(my_units))
    
  } else {
    
    my_units=stakeholder_units
    cols=2
    rows=1
    min_cols=stakeholder_min
    
  }
  
  robustness.plane=SOM_bubbles(SOMlist = SOM.opt, data = robustness.df, units = my_units,
                               title = i, ncol = cols, nrow = rows, my_order=colnames(robustness.df)[-ncol(robustness.df)],
                               min_cols = min_cols, hex = T, text_size = 4.5)
  robustness.planes[[i]]=robustness.plane
  
  if(i == "90% maximin"){ # % is invalide in filename
    i="90maximin"
  }
  if(i == "90% regret from best"){ # % is invalide in filename
    i="90 regret from best"
  }
  
  if(i== "Stakeholder satisficing metrics"){
    
    w=10.6; h=3
    
  }
  
  
  #ggsave(filename = paste0("hexagons - ", i, ".pdf"), plot = robustness.plane, device = "pdf", width = w, height = h )
  
  print(robustness.plane)
  
}


### radar plots (Appendix A11-13 part a)


for(i in names(robustness.list)){
  
  temp.robust=robustness.list[[i]]
  
  if(i=="Stakeholder satisficing metrics"){
    next # cannot plot radar plot of only two metrics
    colnames(temp.robust)=c("Delivery", "Storage")

  } else {
    temp.robust=temp.robust[,-1] # remove policy column
    my_order=c("P.WYR", "LF.Deficit", "LB.Dur", "LB.Max", "LB.Freq","LB.Avg", "M1000", "P3490")
    temp.robust=temp.robust[,my_order]
    
    scale1.robust=data.frame(apply(X=temp.robust, MARGIN = 2, FUN = scale1)) # scale all values 0 to 1
    
    
    
  }
  
  scale1.robust$node=SOM.opt$unit.classif
  radar_grid(x=5, y=3, data = scale1.robust, color = "lightblue", hex_shift = "right", vlcex=1)
  #radar_grid(x=5, y=3, data = scale1.robust, color = "lightblue", hex_shift = F, vlcex=1)
  
  
  mtext(text=i, outer=T, cex=1.3, line = -1.5,adj = 0.415)
  
}

### box plots to take place of satisficing metrics since they are only 2 metrics
### Figure 8, Appendix A14

# need to establish the order to produce the plots
# base R creates plots from top left to bottom right, but the neuron map from Kohonen is labeled from bottom left to top right

n=x*y

node_order=vector() # preallocate
for(r in 1:y){
  add_vec=(n-x*r+1):(n-x*(r-1))
  node_order=c(node_order, add_vec)
  
}

sat.df=robustness.list$`Stakeholder satisficing metrics`
sat.df$node=SOM.opt$unit.classif

box.list=list()
c=1
for (i in node_order){
  
  filter.node=filter(sat.df, node==i)
  
  if(nrow(filter.node)==0){
    temp=NULL
  } else {
    
    filter.node=tidyr::pivot_longer(data = filter.node, cols = c(1,2), names_to="Stakeholder", values_to="Satisficing")
    
    temp=ggplot(data=filter.node, aes(x=Stakeholder, y=Satisficing, fill=Stakeholder))+
      geom_boxplot(outlier.colour = "red", outlier.size = 1)+
      stat_summary(fun=mean, geom="point",shape=18, size=3, color="purple", fill="black")+
      xlab(NULL)+ylab(NULL)+
      ylim(0,1)+
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    temp=temp+ggtitle(i)
  }
  
  if((is.odd(1+floor(c/(x+0.001))) )){ # Shift odd rows to the right
    temp=temp+
      theme(plot.margin = margin(t=0, r=0, b=0, l=3, unit='cm'))
  } else if ((!is.odd(1+floor(i/(x+0.001))) )){ # Shift even rows to the left
    temp=temp+
      theme(plot.margin = margin(t=0, r=3, b=0, l=0, unit='cm'))
  }
  
  if(!(c %in% c(1,6,11))){ # not a far left neuron
    
    temp=temp+
      theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
    
  }
  
  box.list[[c]]=temp
  c=c+1
}


box_plots=ggarrange(plotlist = box.list, nrow=y, ncol=x, common.legend = T, legend="bottom")

box_plots=annotate_figure(box_plots, bottom = text_grob("", vjust=-8.5), left=text_grob("Satisficing fraction", rot=90, y=.6), top= text_grob("Satisficing per neuron", face="bold"))

box_plots


######################################################################################
################# operation diagram topology maps #################################################

# Import Lake Mead DV
Archive='Archive_463_Condensed.txt'
DV=read.table(here("case study data",Archive), header = T, sep = "")

filter.DV=DV

filter.DV$node=SOM.opt$unit.classif

keep=1:463

filter.DV=filter.DV[keep,]

# need to establish the order to produce the plots
# base R creates plots from top left to bottom right, but the neuron map from Kohonen is labeled from bottom left to top right

n=x*y

node_order=vector() # preallocate
for(r in 1:y){
  add_vec=(n-x*r+1):(n-x*(r-1))
  node_order=c(node_order, add_vec)
  
}

## Figure 6

hist.list=list()
c=1
for (i in node_order){
  
  filter.node=filter(filter.DV, node==i)
  
  if(nrow(filter.node)==0){
    temp=NULL
  } else {
    #temp=DV_plot(to_plot = as.numeric(sort(row.names(filter.node))),ID_label = NULL, x=5, y=3, i=i, hex_shift="right")
    temp=DV_plot(to_plot = as.numeric(sort(row.names(filter.node))),ID_label = NULL, x=5, y=3, i=i, hex_shift=F, stat_size = 4, stat_adj = 11)
    
    temp=temp+ggtitle(i)
  }
  
  if(!(c %in% c(1,6,11))){ # not a far left neuron, remove y axis
    
    temp=temp+
      theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
    
  }
  
  
  temp
  
  hist.list[[c]]=temp
  c=c+1
}


combined=ggarrange(plotlist = hist.list, nrow=y, ncol=x, common.legend = T, legend="bottom")

combined=annotate_figure(combined, bottom = text_grob("", vjust=-8.5), left=text_grob("pool elevation [ft msl]", rot=90, y=.6), top= text_grob("Lake Mead decision variables", face="bold"))

combined

#ggsave(filename = "stacked histogram with stats NB.pdf", plot = combined, device = "pdf", width = 11, height = 9 )


### bar plots of loading scores, used in Figures 5 and 8

E1$direction=ifelse(E1$loading>0, "increase", "decrease")

LoadingPC1=ggplot(data=E1, aes(x=forcats::fct_inorder(objective), y=abs(loading), fill=direction, linetype=direction))+
  geom_bar(stat="identity", color="black", size=1)+
  scale_fill_manual(values=c("increase"="mediumorchid2", "decrease"="green4"))+
  xlab("Objective")+
  ylab("loading magnitude")+
  ggtitle("Horizontal axis loading scores")+
  theme(legend.position = c(0.93, .7),
        legend.background = element_rect(fill = "transparent", color = "transparent"), legend.title = element_blank())


LoadingPC1

E2$direction=ifelse(E2$loading>0, "increase", "decrease")

LoadingPC2=ggplot(data=E2, aes(x=forcats::fct_inorder(objective), y=abs(loading), fill=direction, linetype=direction))+
  geom_bar(stat="identity", color="black", size=1)+
  scale_fill_manual(values=c("increase"="mediumorchid2", "decrease"="green4"))+
  xlab("Objective")+
  ylab("loading magnitude")+
  ggtitle("Vertical axis loading scores")+
  theme(axis.title.y=element_text(angle=-90))+ # rotate so legible when plot flipped vertically
  theme(legend.position = c(0.93, .7),
        legend.background = element_rect(fill = "transparent", color = "transparent"), legend.title = element_blank())

LoadingPC2
