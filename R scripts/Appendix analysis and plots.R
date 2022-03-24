# Appendix analysis and plots
# Nathan Bonham
# March 2022

# load packages and custom functions
source(here("R scripts","Library.R")) 

#######################################################################################
###################### scattermatrix of SOW ensemble in Appendix A4 ##############

subsample=readRDS(here("case study data", "SOW ensemble 500_300_100.rds"))

chosen=subsample[[1]]$SOW
scattermat=dplyr::select(chosen, Demand, Mead.PE, Powell.PE, mean)

pm=my_ggpairs(scattermat)

pm=pm + theme(axis.text.x = element_text(angle = 90, hjust = 1))

pm

##################################################################
####### Cumulative Natural Flow above Lees Ferry, AZ #############

# Annual Flow Duration Curves for streamflow time series in the SOW ensemble (this plot takes a minute to render!)
# plot A4, b
plotFDC(subsample = subsample[1], directory = here())

##############################################################################################
##########################  topo and quant error for robustness metrics ######################

# See Appendix A7

### create your own SOM objects with required components to use somQuality

fit.df=data.frame(matrix(ncol=5, nrow=1))
colnames(fit.df)=c("metric", "quant.e", "perc.var", "topo.e", "method" )

robustness.list=readRDS(here('case study data', 'robustness list.rds'))

robustness.list[["maximin"]]=NULL # we removed maximin because worst case SOW less stable

names(robustness.list)=c("mean", "90% maximin", "90% regret", "satisficing")

SOM.opt=readRDS(here('case study data', 'SOMopt.rds'))

x=5; y=3; save_radius=0.540085 # selected params from grid search 
distance.fnc="sumofsquares"; toroidal=F; neighbor.fnc="bubble"

for(i in names(robustness.list)){
  
  mySOMlist=list()
  
  # som grid
  
  mySOMlist[["grid"]]=SOM.opt$grid
  
  # som codes. The codebook vector provides the centroid of each neuron. Here, I will use the value of each objective averaged within each neuron
  
  if(i == "satisficing"){
    
    robust.df=robustness.list[[i]]
    robust.df=scale(robust.df)
    robust.df=data.frame(robust.df)
    robust.df$node=SOM.opt$unit.classif
    avg.df=robust.df %>% group_by(node) %>% summarise(Delivery=mean(Delivery), Storage=mean(Storage))
    
  } else {
    
    robust.df=robustness.list[[i]][,-1] # remove policy column
    robust.df=scale(robust.df)
    robust.df=data.frame(robust.df)
    
    original_order=c("M1000", "LB.Dur", "LB.Freq", "LB.Avg", "LB.Max", "P3490", "P.WYR", "LF.Deficit")
    robust.df=robust.df[,original_order]
    robust.df$node=SOM.opt$unit.classif
    avg.df=robust.df %>% group_by(node) %>% summarise(M1000=mean(M1000), LB.Dur=mean(LB.Dur),
                                                      LB.Freq=mean(LB.Freq), LB.Avg=mean(LB.Avg),
                                                      LB.Max=mean(LB.Max), P3490=mean(P3490),
                                                      P.WYR=mean(P.WYR), LF.Deficit=mean(LF.Deficit))
    
  }
  
  
  codes.list=list()
  codes.list[[1]]=as.matrix(avg.df[,-1])
  
  mySOMlist[["codes"]]=codes.list
  
  # unit classifications
  
  mySOMlist[["unit.classif"]]=SOM.opt$unit.classif
  
  temp=aweSOM::somQuality(mySOMlist, traindat = robust.df[,-ncol(robust.df)])
  optFit=c(i,temp$err.quant, temp$err.varratio, temp$err.topo, "superposition")
  
  fit.df=rbind(fit.df, optFit)
  ###### compare to SOM fit to robustness metrics
  
  X=as.matrix(robust.df[,-ncol(robust.df)])
  temp.eigenVectors=eigen(cov(X))
  RM=temp.eigenVectors$vectors[,1:2] # rotation matrix, which is the eigenvectors corresponding to first and second largest eigenvalues
  PCs=X %*% RM
  PC1range=range(PCs[,1]) # will use these to sample uniformly along PC1 and PC2 within the for loop
  PC2range=range(PCs[,2])
  D1=seq(PC1range[1], PC1range[2], length.out = x) # sequence of nuerons along PC1
  D2=seq(PC2range[1], PC2range[2], length.out = y) # sequence of nuerons along PC1
  
  PG=expand.grid(D1,D2) # create rectangular matrix where D1 is repeated for every value of D2, projected in PC space
  
  # plot(x=PG$Var1,y=PG$Var2, type = "p" ) # to see example of the grid in PC space
  # points(PCs[,1], PCs[,2], col="red") # to see the policies in PC space
  # unproject the projected grid (PG) back into original data space
  
  IG=as.matrix(PG) %*% t(RM) # neuron initialization matrix
  
  
  radius=fraction2radius(fraction = save_radius , x = x, y = y, shape = "hexagonal")
  
  temp.som=som(X =X,radius= radius, dist.fcts=distance.fnc,
               grid=somgrid(xdim=x, ydim=y, topo = "hexagonal", toroidal = toroidal, neighbourhood.fct = neighbor.fnc), rlen=10000, keep.data=T,
               mode="batch", init=IG)
  
  temp=aweSOM::somQuality(temp.som, traindat = X)
  SOMfit=c(i,temp$err.quant, temp$err.varratio, temp$err.topo, "SOM fit")
  fit.df=rbind(fit.df, SOMfit)
  
}

#fit.df=rbind(fit.df, firstSOMfit) # add metrics of SOM fit to optimization objectives
fit.df=fit.df[-1,] # remove initialization row

fit.df[,2]=as.numeric(fit.df[,2])
fit.df[,3]=as.numeric(fit.df[,3])
fit.df[,4]=as.numeric(fit.df[,4])
fit.df$topo.correct=1-fit.df$topo.e

fit.df$frac.var=fit.df$perc.var/100

fit.df=fit.df[fit.df$metric!="satisficing",] # remove satisficing from bar plot

### plot, Appendix A10

size=1.2; lt=3
my_angle=15
hjust=.5
vjust=.9

quant=ggplot(data=fit.df, aes(x=fct_inorder(metric), y=quant.e, fill=method))+
  geom_bar(stat = "identity", color="black", position=position_dodge())+
  coord_cartesian(ylim = c(0,2))+
  ggtitle("Quantization error")+
  ylab("error")+
  xlab("metric")+
  theme(axis.text.x = element_text(angle = my_angle, hjust=hjust, vjust=vjust))+
  scale_linetype_manual(name="", values = c(lt), 
                        guide = guide_legend(override.aes = list(color = c("black"), size=size)))

perc.var=ggplot(data=fit.df, aes(x=fct_inorder(metric), y=perc.var, fill=method))+
  geom_bar(stat = "identity", color="black", position=position_dodge())+
  coord_cartesian(ylim = c(70,100))+
  ggtitle("Percent of variance explained")+
  ylab("percent")+
  xlab("metric")+
  theme(axis.text.x = element_text(angle = my_angle, hjust=hjust, vjust=vjust))+
  scale_linetype_manual(name="", values = c(lt), 
                        guide = guide_legend(override.aes = list(color = c("black"), size=size)))

frac.var=ggplot(data=fit.df, aes(x=fct_inorder(metric), y=frac.var, fill=method))+
  geom_bar(stat = "identity", color="black", position=position_dodge())+
  coord_cartesian(ylim = c(0,1))+
  ggtitle("Fraction of variance explained")+
  ylab("fraction")+
  xlab("metric")+
  theme(axis.text.x = element_text(angle = my_angle, hjust=hjust, vjust=vjust))+
  scale_linetype_manual(name="", values = c(lt), 
                        guide = guide_legend(override.aes = list(color = c("black"), size=size)))



topo.e=ggplot(data=fit.df, aes(x=fct_inorder(metric), y=topo.e, fill=method))+
  geom_bar(stat = "identity", color="black", position=position_dodge())+
  ylim(0,1)+
  ggtitle("Topographic error")+
  ylab("error")+
  xlab("metric")+
  theme(axis.text.x = element_text(angle = my_angle, hjust=hjust, vjust=vjust))+
  scale_linetype_manual(name="", values = c(lt), 
                        guide = guide_legend(override.aes = list(color = c("black"), size=size)))

topo.correct=ggplot(data=fit.df, aes(x=fct_inorder(metric), y=topo.correct, fill=method))+
  geom_bar(stat = "identity", color="black", position=position_dodge())+
  coord_cartesian(ylim = c(0,1))+
  ggtitle("Topographic skill: 1 - topo error")+
  ylab("")+
  xlab("metric")+
  theme(axis.text.x = element_text(angle = my_angle, hjust=hjust, vjust=vjust))+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  scale_linetype_manual(name="", values = c(lt), 
                        guide = guide_legend(override.aes = list(color = c("black"), size=size)))

ggarrange(frac.var, topo.correct, nrow=1, common.legend = T, legend = "right" )

