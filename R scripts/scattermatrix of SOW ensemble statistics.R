# Plot scatter matrix of SOW ensemble for post-MORDM paper
# Nathan Bonham
# January 2022
library(dplyr)
library(GGally)

subsample=readRDS("G:/My Drive/CU Boulder/CRB publications/Uncertainty characterization and RDM sensitivity paper/R/data/SOW ensemble 500_300_100.rds")

chosen=subsample[[1]]$SOW
scattermat=dplyr::select(chosen, Demand, Mead.PE, Powell.PE, mean)

my_fn <- function(data, mapping, ...){
  meanX <- mean(eval_data_col(data, mapping$x), na.rm = TRUE)
  
  text=as.character(round(meanX, digits=2))
  #text=paste('mean=',round(meanX, digits=2))
  
  p <- ggplot(data = data, mapping = mapping) + 
    geom_histogram(bins = 20)+
    #geom_density(adjust=2)+
    geom_vline(xintercept = meanX,color='red',linetype='dashed' ) +
    
    annotate('text', label=text,x=Inf, y=Inf, vjust=1, hjust=1, size=4, col='red')
  p
}

pm<-ggpairs(
  scattermat,
  upper=list(continuous="blank"),
  lower=list(continuous=wrap('points', alpha=0.3)),
  diag=list(continuous=my_fn)
  
)
pm=pm + theme(axis.text.x = element_text(angle = 90, hjust = 1))

pm

path="G:/My Drive/CU Boulder/CRB publications/Uncertainty characterization and RDM sensitivity paper/R/Figures"
filename='cLHS N500 scattermate.pdf'
ggsave(filename = filename, plot=pm, device='pdf', path = path, width=5.5, height = 5.5, units = 'in')

