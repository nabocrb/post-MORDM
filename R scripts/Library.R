# Library for post-MORDM
# March 2022
# Nathan Bonham

rm(list=ls()) # clear the environment

# load libraries

library(here) # for relative file paths
library(magick)
library(stringr)
library(dplyr)
library(ggplot2)
library(GGally)
library(plotly)
library(yarrr) #transparent function
library(forcats) # fct_inorder
library(cowplot) # plot_grid
library(aweSOM)
library(ggpubr) # ggarrange
library(ggforce) # circles
library(GGally) # scatter matri plots
library(fmsb) # radar plots
library(shadowtext)
library(kohonen) # SOM
library(aweSOM) # quality metrics: quantization error (perc. variance), topo error
library(lhs) # latin hypercube design of SOM parameters
library(clusterSim) # Davies-Bouldin index

####################### load objectives data from MOEA #####################################

tradeoff_dfs=readRDS(here("case study data", "tradeoff_dataframes.rds"))

optimization=tradeoff_dfs$optimization

######################################### Load SOW ensemble ###################################

temp=readRDS(here("case study data", "SOW ensemble 500_300_100.rds"))
cLHS500all=temp[[1]]$SOW

### select 16 observable uncertainty metrics from cLHS500

cLHS500=dplyr::select(cLHS500all, contains(c("Scenario", "TraceNumber","Demand", ".PE","yr","MaxDur.ter", "MaxDur.hist")))
cLHS500=dplyr::select(cLHS500, -contains("RolDef"))

############################################################################################
####################### find policy closest to prototype vector of each SOM nueron (I call them captains) ###################################################

findCaptains=function(mySOM){ # function not used in post-MORDM paper
  
  ### find policy closest to each node
  
  temp.df=data.frame(node=mySOM$unit.classif, distance=mySOM$distances, policy=1:length(mySOM$distances))
  
  nodeID=vector()
  policyID=vector()
  
  c=1
  for (i in sort(unique(mySOM$unit.classif))){
    
    filter.node=dplyr::filter(temp.df, node==i)
    
    nodeID[c]=i
    policyID[c]=filter.node$policy[which.min(abs(filter.node$distance))]
    
    c=c+1
    
  }
  
  captains=data.frame(node=nodeID, policy=policyID)
  
  return(captains)
}


##############################################################################################
####################### robustness metrics ###################################################

# objective values for each policy and each SOW
obj_all=read.table(here("case study data", "objectives_all463.txt")) # load data for calculating robustness

# regret2 is short for regret from best

regret2=function(data=obj, objectives=c('LB.Shortage.Volume', 'Mead.1000', 'Powell.3490'),
                 best_if=c('min', 'min', 'min'), 
                 policy_ID_column=2, SOW_column=1, SOW_agg_method='percentile',
                 percentile=0.9, Obj_agg_method='none', scale=T){
  
  
  ############ Calculate best and worst performance in every SOW, for each objective
  
  SOW_iter=unique(data[,SOW_column])
  regret2_data=data[,c(policy_ID_column,SOW_column,which(colnames(data) %in% objectives))] # remove objectives we are not using
  
  best=matrix(nrow=length(SOW_iter), ncol=ncol(regret2_data)-2) # preallocated matrix to store best performance values of objective i in SOW j
  worst=matrix(nrow=length(SOW_iter), ncol=ncol(regret2_data)-2) # preallocated matrix to store worst performance values
  
  # find optimal performance for objective i in SOW j. Also find worst performance for objective i in SOW j to calculate max deviation later  
  ######### You should remove this part from the policy loop. Redundant computation #################
  for (i in 1:length(SOW_iter)){
    filter.SOW=regret2_data[which(regret2_data[,2]==SOW_iter[i]),]
    
    for (c in 1:ncol(best)){ # loop through objectives, find min or max as best performance
      
      if (best_if[c]=='min'){
        best[i,c]=min(filter.SOW[,(2+c)])
        worst[i,c]=max(filter.SOW[,(2+c)])
      } else if (best_if[c]=='max'){
        best[i,c]=max(filter.SOW[,(2+c)])
        worst[i,c]=min(filter.SOW[,(2+c)])
      } else {
        paste('ERROR: this function only supports min or max for finding best objective performance in regret type II calculations')
      }
      
    } # end objectives loop
    
    
  } # end SOW iteration
  
  
  ############## end of computing best and worst matrices ###############
  
  ####### prepare data frame to store results
  
  if (Obj_agg_method=='none'){
    ncol=length(objectives)+1 # one per objective plus one for policy ID
  } else { # use default aggregation (normalize and sum)
    ncol=2 # one column for satisficing plus one for policy ID
  }
  
  nrow=length(unique(data[, policy_ID_column])) # one row per policy
  
  metric.df=matrix(NA, nrow=nrow, ncol=ncol)
  metric.df=data.frame(metric.df)
  colnames(metric.df)=c('policy', 'regret2')
  metric.df[,1]=unique(data[, policy_ID_column])
  
  policy_iter=unique(data[, policy_ID_column]) # id of policies to loop through
  
  
  for (r in 1:length(policy_iter)){ ########## loop through policies
    
    filter.policy=data[which(data[,policy_ID_column]==policy_iter[r]),]
    
    ##################### Compute type 2 regret #################################
    
    regret2_calcs=matrix(ncol=length(objectives)+1, nrow=nrow(filter.policy))
    
    
    for (cr in 1:length(objectives)){ # loop through criteria
      
      numerator=abs(filter.policy[,which(colnames(filter.policy)==objectives[cr])] - best[,cr])
      # denominator=filter.policy[,which(colnames(filter.policy)==objectives[cr])] # I don't like this as the normalizing factor
      
      ##### max deviation across all policies for each SOW #####################
      
      max_deviation=best[,cr] - worst[,cr] # this is the maximum deviation across all POLICIES for a given SOW (the rows). Ie, best performance in SOW j across all policies minus worst performance in SOW j across all policies
      
      # max_deviation=range(data[,which(colnames(data)==objectives[cr])])[2]-range(data[,which(colnames(data)==objectives[cr])])[1] # this isn't SOW specific
      denominator=abs(max_deviation) # scale everything zero to 1 by dividing by the largest deviation possible across all Policy and SOW
      
      if(scale==F){denominator=1}
      
      regret2_calcs[,cr]=numerator/denominator
      
      denom_0=which(round (denominator, 2)==0) # anywhere the demoninator is zero, the result is NaN.
      
      regret2_calcs[denom_0,cr]=numerator[denom_0] # if denominator is zero, just use numerator
      
      
    } # end regret2 criteria loop  
    
    if (length(objectives)==1){
      regret2_calcs[, length(objectives)+1]=regret2_calcs[,1]
     } else {
      regret2_calcs[, length(objectives)+1]=rowSums(regret2_calcs[,1:length(objectives)])
    }
    
    
    ######### SOW aggregation ################
    
    if (SOW_agg_method=='percentile'){
      
      metric.df$regret2[r]=quantile(regret2_calcs[,ncol(regret2_calcs)],percentile)
      
    } else {
      agg_func=match.fun(SOW_agg_method) # transform string into function. 
      metric.df$regret2[r]=agg_func(regret2_calcs[,ncol(regret2_calcs)])
      
    }
    
    ######## end SOW aggregation #############
    
  } # end policy loop
  
  best=data.frame(best)
  worst=data.frame(worst)
  colnames(best)=objectives
  colnames(worst)=objectives
  
  return(list(regret2=metric.df, best_per_SOW=best, worst_per_SOW=worst))
  
} # end function


######################################## satisficing #####################################

satisficing=function(data=obj, objectives=c('Lee.Ferry.Deficit', 'Mead.1000', 'Powell.3490'),
                     thresholds=c(0, 10, 1), fail_if_inequality=rep('greater', length(objectives)), 
                     n_satisficing=length(objectives), policy_ID_column=2, SOW_column=1){
  
  ####### prepare data frame to store results
  
  ncol=2 # one column for satisficing plus one for policy ID
  nrow=length(unique(data[, policy_ID_column])) # one row per policy
  
  metric.df=matrix(NA, nrow=nrow, ncol=ncol)
  metric.df=data.frame(metric.df)
  colnames(metric.df)=c('ID', 'satisficing')
  metric.df[,1]=unique(data[, policy_ID_column])
  
  policy_iter=unique(data[, policy_ID_column]) # id of policies to loop through
  
  if(n_satisficing > length(objectives)){
    
    n_satisficing=length(objectives)
    print('n_satisficing must be <= number of objectives. n_satisficing has been changed to length(objectives).')
    
  }
  
  if(length(fail_if_inequality) != length(objectives)){
    
    fail_if_inequality=rep('greater', length(objectives))
    print('length(fail_if_inequality) must = length(objectives). fail_if_inequality has been changed to rep(greater, length(objectives)).')
    
  }
  
  
  for (r in 1:length(policy_iter)){ ########## loop through policies
    
    filter.policy=data[which(data[,policy_ID_column]==policy_iter[r]),]
    
    ##################### Compute satisficing #################################
    
    satisficing_calcs=matrix(ncol=length(objectives)+2, nrow=nrow(filter.policy))
    
    for (cr in 1:length(objectives)){
      
      if (fail_if_inequality[cr]=='greater'){
        
        satisficing_calcs[,cr]=ifelse(filter.policy[,which(colnames(filter.policy)==objectives[cr])] > thresholds[cr], 0,1) # a 1 means 'met criteria'
        
      } else if (fail_if_inequality[cr]=='less'){
        
        satisficing_calcs[,cr]=ifelse(filter.policy[,which(colnames(filter.policy)==objectives[cr])] < thresholds[cr], 0,1) # a 1 means 'met criteria'
        
      } else {
        
        satisficing_calcs[,cr]=NA
        
      }
      
      
    } # satisficing criteria loop 
    
    if (length(objectives)==1){ # only one objective, don't need to do rowSums operator
      satisficing_calcs[, length(objectives)+1]=satisficing_calcs[,1]
      satisficing_calcs[, length(objectives)+2]=ifelse(satisficing_calcs[, length(objectives)+1] < n_satisficing , 0, 1) # a 1 in final column means yes, satisficing
    } else { # more than one objective. Use rowSums
      satisficing_calcs[, length(objectives)+1]=rowSums(satisficing_calcs[,1:length(objectives)])
      satisficing_calcs[, length(objectives)+2]=ifelse(satisficing_calcs[, length(objectives)+1] < n_satisficing , 0, 1) # a 1 in final column means yes, satisficing
    }
    
    metric.df$satisficing[r]=sum(satisficing_calcs[,ncol(satisficing_calcs)])/nrow(satisficing_calcs)
    
    
  } # end policy loop
  
  return(metric.df)
  
} # end function


# default percentiles assume minimization. If maximization, change percentile to 0

maximin=function(data=obj, objectives=c('LB.Shortage.Volume', 'Mead.1000', 'Powell.3490'),
                 percentile=rep(100, length(objectives)), 
                 policy_ID_column=2, SOW_column=1){
  
  ####### prepare data frame to store results
  
  ncol=length(objectives)+1 #plus one for policy ID
  nrow=length(unique(data[, policy_ID_column])) # one row per policy
  
  metric.df=matrix(NA, nrow=nrow, ncol=ncol)
  metric.df=data.frame(metric.df)
  colnames(metric.df)=c('policy', objectives)
  metric.df[,1]=unique(data[, policy_ID_column])
  
  policy_iter=unique(data[, policy_ID_column]) # id of policies to loop through
  
  
  for (r in 1:length(policy_iter)){ ########## loop through policies
    
    filter.policy=data[which(data[,policy_ID_column]==policy_iter[r]),]
    
    ##################### Compute maximin #################################
    
    
    calcs=matrix(ncol=length(objectives), nrow=nrow(filter.policy))
    
    for (cr in 1:length(objectives)){
      
      
      metric.df[r,(cr+1)]=quantile( filter.policy[,which(colnames(filter.policy)==objectives[cr])], percentile[cr]/100) #return the desired quantile. For case of minimization and maximin, want 100%
      
    } #  criteria loop  
    
    
    
  } # end policy loop
  
  return(metric.df)
  
} # end function


LaplacePIR=function(data=obj, objectives=c('LB.Shortage.Volume', 'Mead.1000', 'Powell.3490'),
                    policy_ID_column=2, SOW_column=1){
  
  ####### prepare data frame to store results
  
  ncol=length(objectives)+1 #plus one for policy ID
  nrow=length(unique(data[, policy_ID_column])) # one row per policy
  
  metric.df=matrix(NA, nrow=nrow, ncol=ncol)
  metric.df=data.frame(metric.df)
  colnames(metric.df)=c('policy', objectives)
  metric.df[,1]=unique(data[, policy_ID_column])
  
  policy_iter=unique(data[, policy_ID_column]) # id of policies to loop through
  
  
  for (r in 1:length(policy_iter)){ ########## loop through policies
    
    filter.policy=data[which(data[,policy_ID_column]==policy_iter[r]),]
    
    ##################### Compute metric #################################
    
    
    for (cr in 1:length(objectives)){
      
      metric.df[r, (cr+1)]=mean(filter.policy[,which(colnames(filter.policy)==objectives[cr])])
      
    } #  objective loop  
    
    
    
  } # end policy loop
  
  return(metric.df)
  
} # end function

##################################################################
############### SOM helper functions #############################

#### convert fraction to SOM radius

fraction2radius=function(fraction, x,y,shape){
  
  # I took this code from the source code of the supersom function in Kohonen package and turned it into a function
  # to calculate the radius given a fraction of the max node-to-node distance
  
  grid = somgrid(x,y,shape)
  grid <- check.somgrid(grid)
  nhbrdist <- unit.distances(grid)
  radius = quantile(nhbrdist, fraction)
  
  return(radius)
  
}

check.somgrid=function (grd) # taken from Kohonen source code
{
  mywarn <- FALSE
  if (is.null(grd$toroidal)) {
    mywarn <- TRUE
    grd$toroidal <- FALSE
  }
  if (is.null(grd$neighbourhood.fct)) {
    mywarn <- TRUE
    grd$neighbourhood.fct <- factor("bubble", levels = c("bubble", 
                                                         "gaussian"))
  }
  if (mywarn) 
    warning("Added defaults for somgrid object - ", "you are probably using the somgrid function ", 
            "from the class library...")
  grd
}

unit.distances=function (grid, toroidal) # taken from Kohonen source code
{
  if (missing(toroidal)) 
    toroidal <- grid$toroidal
  if (!toroidal) {
    if (grid$topo == "hexagonal") {
      return(as.matrix(stats::dist(grid$pts)))
    }
    else {
      return(as.matrix(stats::dist(grid$pts, method = "maximum")))
    }
  }
  np <- nrow(grid$pts)
  maxdiffx <- grid$xdim/2
  maxdiffy <- max(grid$pts[, 2])/2
  result <- matrix(0, np, np)
  for (i in 1:(np - 1)) {
    for (j in (i + 1):np) {
      diffs <- abs(grid$pts[j, ] - grid$pts[i, ])
      if (diffs[1] > maxdiffx) 
        diffs[1] <- 2 * maxdiffx - diffs[1]
      if (diffs[2] > maxdiffy) 
        diffs[2] <- 2 * maxdiffy - diffs[2]
      if (grid$topo == "hexagonal") {
        result[i, j] <- sum(diffs^2)
      }
      else {
        result[i, j] <- max(diffs)
      }
    }
  }
  if (grid$topo == "hexagonal") {
    sqrt(result + t(result))
  }
  else {
    result + t(result)
  }
}

###########################################################################################
######################################### Plotting ########################################

#### topology map plotting functions

### radar plots 

scale1=function(data){
  max=max(data)
  min=min(data)
  
  temp=(data-min)/(max-min)
  
  return(ifelse(temp>1, 1, temp))
  
}


# https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
my_radar<- function(data, color = "#00AFBB", alpha=0.5, line_color="black",
                    vlabels = colnames(data), vlcex = 0.7, plty=1,
                    caxislabels = c(0, .25, .5, .75, 1), title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol =  scales::alpha(line_color, alpha), pfcol = scales::alpha(color, alpha), plwd = 2, plty = plty,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.3,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = NULL, ...
  )
  
  mtext(text = title, line = -1.5, adj=.3, font = 2)
  
}

is.odd=function(x){
  
  return(x%%2 != 0)
  
}

radar_grid=function(x,y,data, color, alpha=0.3, line_color="black", hex_shift="right", vlcex=1.3, vlabels=colnames(data), caxislabels = c(0, .25, .5, .75, 1), title=T){
  
  par(mfrow = c(y,x))
  par(mar = c(.05, 2, .05, 2), oma=c(0,0,0,0))
  par(mar = c(.00, 0, 0, 0), oma=c(0,0,0,0))
  
  
  my_line="solid"
  
  # need to establish the order to produce the plots
  # base R creates plots from top left to bottom right, but the neuron map from Kohonen is labeled from bottom left to top right
  
  n=x*y
  
  node_order=vector() # preallocate
  for(r in 1:y){
    add_vec=(n-x*r+1):(n-x*(r-1))
    node_order=c(node_order, add_vec)
    
  }
  
  for(i in node_order){
    
    
    
    filter.node=filter(data, node==i)
    
    alpha=0.3
    
    if(nrow(filter.node)==0){
      
      filter.node[1,]=999
      alpha=0
      # my_line="blank"
    }
    
    # data must have max and min values for each variable in row 1 and 2 for fmsb radar chart. See here:
    # https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
    
    range=rbind(rep(1, ncol(data)), rep(0, ncol(data)))
    
    colnames(range)=colnames(data)
    range=data.frame(range)
    range$node=c(999, 999) # 999 for max and min rows
    
    temp=rbind(range, filter.node)[,-ncol(data)] # add the max and min rows, remove node column
    
    if(hex_shift=="right" & (is.odd(1+floor(i/(x+0.001))) )){ # Shift odd rows to the right
      par(mar=c(.05,6, .05, 0)) # bottom, left, top, right. I use x+0.001 for the case where i is a multiple of x 
    } else if (hex_shift=="right" & (!is.odd(1+floor(i/(x+0.001))) )){ # Shift even rows to the left
      par(mar=c(.05,0, .05, 6)) # bottom, left, top, right
    }
    
    if(title){my_title=i} else { my_title = ""}
    
    my_radar(data=temp, color =rep("lightblue", nrow(filter.node)), alpha = alpha, plty=my_line, title = my_title, vlcex=vlcex, vlabels=vlabels, caxislabels = caxislabels)
    
    
    
  }
  
  
}

### Hexagonal or rectangular grid SOM topology maps 

# my method for defining the coordinates of the SOM grids taken from: http://blog.schochastics.net/post/soms-and-ggplot/

SOM_bubbles=function(SOMlist, data, units, my_order, title="", ncol, nrow, min_cols=rep(T, ncol(data)-1), hex=F, hex_size=0.05, no_fill=F, labels=T, decimals=2, text_size=2.5){
  
  nodeAVG= data %>% group_by(node) %>% summarise_each(list(mean=mean))
  coords=SOMlist$grid$pts # http://blog.schochastics.net/post/soms-and-ggplot/
  
  if (nrow(nodeAVG)!=nrow(coords)){ # happens when a node has been eliminated by non-domination filter
    
    all_nodes=1:nrow(coords)
    kept_nodes=sort(unique(nodeAVG$node))
    
    missing_node=all_nodes[which(!(all_nodes %in% kept_nodes))]
    
    for(r in missing_node){
      add_row=c(r, rep(NA, ncol(nodeAVG)-1))
      nodeAVG=rbind(nodeAVG, add_row)
      
    }
    
    nodeAVG=dplyr::arrange(nodeAVG, node)
    
  }
  
  nodeAVG$x=coords[,1]
  nodeAVG$y=coords[,2]  
  
  circle_plot.list=list()
  
  for(i in 1:length(my_order)){
    
    label=my_order[i]
    my_var=colnames(nodeAVG)[i+1] # plus one to skip node column
    my_unit=units[i]
    
    
    if(min_cols[i]){ # if a minimization column, blue to brown gradient
      #my_cols=rev(c("brown", "Lightgreen", "DeepSkyBlue"))
      my_cols=c("gray10", "white")
    } else { # max column, use brown to blue gradient
      #my_cols=c("brown", "Lightgreen", "DeepSkyBlue")
      my_cols=rev(c("gray10", "white"))
    }
    
    if(units[i]=="fraction SOW"){
      my_range=c(min(dplyr::select(data, -node)),max(dplyr::select(data, -node)))
    } else {
      my_range=NULL
    }
    
    
    if(labels){
      
      nodeAVG[[my_var]]=round(nodeAVG[[my_var]], decimals)
      #nodeAVG[["negative"]]=-1*nodeAVG[[my_var]]
      
    }
    
    
    my_plot=ggplot(data=nodeAVG)
    
    if(hex==F){
      my_plot=my_plot+
        geom_circle(aes_string(x0="x", y0="y", r=0.4, fill=my_var))+
        #geom_circle(aes_string(x0="x", y0="y", r=0.4), fill="transparent")+ # to create colorless map
        coord_fixed()+
        theme(plot.margin=margin(t=-100, r=5, b=-100, l=5))
    } else {
      my_plot=my_plot+
        stat_summary_hex(aes_string(x="x", y="y", z=my_var), binwidth = hex_size, color="black")+
        theme(plot.margin=margin(t=25,r=0,b=0,l=0))+
        coord_cartesian(xlim=c(min(nodeAVG$x)-.5,max(nodeAVG$x)+.5), ylim=c(min(nodeAVG$y)-.5, max(nodeAVG$y)+.5))
      
      # nice example of hexagon plot https://andrewpwheeler.com/2019/08/07/making-a-hexbin-map-in-ggplot/
    }
    
    if(no_fill){
      my_cols=c("white", "white")
    }
    
    if(labels){
      
      my_plot=my_plot+
        geom_shadowtext(aes_string(x="x", y="y", label=my_var), color="white", size=text_size) #, color="dodgerblue1"
    }
    
    my_plot=my_plot+  
      #scale_fill_continuous(limits=range(nodeAVG[[my_var]]), low="blue", high="red", name=my_unit)+
      scale_fill_gradientn(colors=my_cols, name=my_unit, limits=my_range)+
      #scale_color_gradientn(colors=rev(my_cols), limits=my_range)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
      ggtitle(label)+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_blank() )
    
    circle_plot.list[[label]]=my_plot
    
  }
  
  
  circles_per_obj=ggarrange(plotlist = circle_plot.list, ncol = ncol, nrow=nrow )
  
  circles_per_obj=annotate_figure(circles_per_obj, top= text_grob(title, face="bold", vjust=2, size = 16))
  
  return(circles_per_obj)
  
}

### Reservoir operation diagrams

# import Lake Mead DV values

bar_plot_data=readRDS(here("case study data",'data for stacked bar plot.rds'))
long_data=bar_plot_data$long_data
wide_data=bar_plot_data$wide_data
wide_data$ID=1:nrow(wide_data)


# function to change ggplot legend size: https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend <- function(myPlot=bar_plot, pointSize = 0.25, textSize = 8, spaceLegend = 0.03) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_blank(), #element_text(size = textSize)
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

DV_plot=function(long.data=long_data, wide.data=wide_data, metric= NULL,
                 to_plot=1:10, xlabel='', ylabel='', ID_label="ID", preferred_direction='max', y_axis2=F,
                 labelsize=3, shrink_legend=F, interactive=F, volume_labs=F, v_lab_nudge=-4, hex_shift=NA, x=5, y=3, i=i, summary_stats=T, stat_size=3,stat_adj=10){
  
  # filter for chosen policies
  filter.long=dplyr::filter(long.data, policy %in% to_plot)
  filter.wide=dplyr::filter(wide.data, ID %in% to_plot)
  
  # get rank, append to data frames
  
  correction=ifelse(preferred_direction == 'min', 1, -1) # to handle metrics that should be minimized and maximizied accordingly in the ranking
  decreasing=ifelse(preferred_direction == 'min', F, T) # needed for add_lines in plot function
  
  if (!is.null(metric)){
    
    filter.metric=dplyr::filter(metric, ID %in% to_plot)
    filter.metric$rank=rank(correction*filter.metric[[metric_label]], ties.method = 'first')
    
  }
  
  filter.wide$rank=if(is.null(metric)){1:length(to_plot)} else {filter.metric$rank}
  filter.long$rank=NA
  for(i in 1:length(filter.long$rank)){ # for loop to avoid errors when different number of tiers
    filter.long$rank[i]=filter.wide$rank[which(filter.wide$ID==filter.long$policy[i])]
  }
  

  ############################# plotting ################################
  text_size=labelsize
  
  n_policies=nrow(filter.wide)
  
  filter.long$volume[filter.long$volume==0]=NA
  
  breaks=seq(0,2400, length.out=13)
  my_labs= c("[0-200)","[200-400)", "[400-600)", "[600-800)", "[800-1000)", "[1000-1200)","[1200-1400)",
             "[1400-1600)","[1600-1800)", "[1800-2000)", "[2000-2200)", "[2200-2400)")
  filter.long$vol_group=cut(filter.long$volume, breaks = breaks, include.lowest = T, right = F, labels = my_labs)
  
  color_func=colorRampPalette(colors=c("greenyellow","lightgoldenrod1", "goldenrod1", "darkorange1", "maroon2","firebrick3"))
  my_cols=color_func(12)
  names(my_cols)=my_labs
  
  
  bar_plot=ggplot()+
    geom_bar(data=filter.long, aes(fill=vol_group,x=rank, y=delta), position = 'stack', stat = 'identity', color='darkgrey')+ 
    # geom_text(data=filter.wide, aes(x=rank, y= policy_lab_y, label=SOM_node), nudge_y = 10, size=text_size, check_overlap = T)+ I have removed SOM node for NOW
    scale_fill_manual(values=my_cols, na.value="black", name="Volume [KAF]", na.translate=F, drop=F)+
    xlab(xlabel)+
    ylab(ylabel)+
    theme(plot.title = element_text(size=10), plot.margin=margin(t=0, r=0, b=0, l=0, unit='pt'))+
    coord_cartesian(ylim=c(905,1110), xlim=c(0.5, min((n_policies+.5),20)))
  
  if(!is.null(ID_label)){
    
    bar_plot=bar_plot+
      geom_text(data=filter.wide, aes_string(x="rank", y= "policy_lab_y", label=ID_label), nudge_y = 4, size=text_size, check_overlap = T)
    
  }
  
  if(!is.na(hex_shift)){
  if(hex_shift=="right" & (is.odd(1+floor(i/(x+0.001))) )){ # Shift odd rows to the right
    bar_plot=bar_plot+
      theme(plot.margin = margin(t=0, r=0, b=0, l=3, unit='cm'))
  } else if (hex_shift=="right" & (!is.odd(1+floor(i/(x+0.001))) )){ # Shift even rows to the left
    bar_plot=bar_plot+
      theme(plot.margin = margin(t=0, r=3, b=0, l=0, unit='cm'))
  }
  }
  
  if(volume_labs){bar_plot=bar_plot+geom_text(data=filter.long, aes(x=rank, y= elevation, label=v_lab),color='black', nudge_y = v_lab_nudge, size=text_size, check_overlap = T)}
  
  if(is.null(metric)){bar_plot=bar_plot+theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())}
  
  if(shrink_legend){bar_plot=addSmallLegend(bar_plot)}
  
  if(summary_stats){
    bar_plot=bar_plot+
      geom_hline(yintercept=mean(filter.wide$T1e), linetype="dashed", color="grey45")+
      annotate("text", x=min((n_policies+.5),20+.5), y=mean(filter.wide$T1e)+2,label=paste0(round(mean(filter.wide$T1e))), vjust=0,hjust=1, size=stat_size)+
      annotate("text", x=min(filter.long$rank)-.5, y=1110,label=paste0(round(mean(filter.wide$T1V)),",", round(mean(filter.wide$maxVol))), hjust=0, size=stat_size)
    #annotate("text", x=min(filter.long$rank)-.5, y=1110,label=paste0("T1e=", round(mean(filter.wide$T1e))), hjust=0, size=stat_size)
    #annotate("text", x=min(filter.long$rank)-.5, y=1110-stat_adj,label=paste0("T1V=", round(mean(filter.wide$T1V))), hjust=0, size=stat_size)+
    #annotate("text", x=min(filter.long$rank)-.5, y=1110-2*stat_adj, label=paste0("maxVol=", round(mean(filter.wide$maxVol))), hjust=0, size=stat_size)#+
    #annotate("text", x=min(filter.long$rank)-.5, y=1080, label=paste0("ntiers=", round(mean(filter.wide$nTiers),1)), hjust=0, size=stat_size)
    
  }
  
  # ggtitle(paste(metric_label, 'rank for selected policies', sep=' '))+
  
  
  if(interactive){
    
    int_plot=ggplotly(p=bar_plot, tooltip = c('rank','elevation', 'volume'), dynamicTicks=T, originalData=F) # convert ggplot to interactive plotly html
    # add legend title in correct location
    int_plot=int_plot %>%   layout(legend = list(
      orientation = "v", title=list(text=" Tier "))
    )
    
    if (y_axis2==T & !is.null(metric)){ # only add second y axis if y_axis2==T and a metric was given
      
      int_plot_2y=int_plot %>%
        add_lines(data=filter.metric, x=~sort(filter.metric$rank), y=~sort(filter.metric[[metric_label]], decreasing = decreasing), yaxis='y2',
                  inherit=FALSE, showlegend=FALSE, line=list(color='purple', width=2, dash='dash')) %>%
        layout(yaxis2 = list(overlaying = "y", side = "right",
                             tickfont = list(color = 'purple', size=10), color = 'purple',
                             title = metric_label),
               legend = list(x = 1.05, y = 0.95), xaxis=list(range=c(0, min((n_policies+1),20))), yaxis=list(range=c(885,1110))
        )
      
    } else { # do not add second axis.
      
      int_plot_2y=int_plot %>%
        layout(legend = list(x = 1.05, y = 0.95), xaxis=list(range=c(0, min((n_policies+1),20)), title=""), yaxis=list(range=c(885,1110)))
      
    }
    
    
    int_plot_2y$x$layout$xaxis$autorange = FALSE # need to tell plotly to NOT change the axis range to fit all data
    int_plot_2y$x$layout$yaxis$autorange = FALSE
    return(int_plot_2y)
    
  } else {
    return(bar_plot) # not interactive
  }
  
  
}

#### additional plotting functions, used for figures in appendix

### ggpairs helper functions (Appendix A4 a)
## function to add vertical line to pdf to indicate mean of each metric plot on diagonal of ggpairs functions
ggpairsMean <- function(data, mapping, ...){ 
  meanX <- mean(eval_data_col(data, mapping$x), na.rm = TRUE)
  
  text=as.character(round(meanX, digits=2))
  #text=paste('mean=',round(meanX, digits=2))
  
  p <- ggplot(data = data, mapping = mapping) + 
    geom_density(adjust=3, kernel="gaussian")+
    geom_vline(xintercept = meanX,color='red',linetype='dashed' ) +
    
    annotate('text', label=text,x=Inf, y=Inf, vjust=1, hjust=1, size=4, col='red')
  p
}

## histograms, instead of pdf

ggpairs_hist <- function(data, mapping, ...){
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

### ggpairs scatter matrix

my_ggpairs=function(chosen){

  pm<-ggpairs(
    chosen,
    upper=list(continuous='cor'),
    lower=list(continuous=wrap('points', alpha=0.3)),
    diag=list(continuous=ggpairs_hist)
    
  )
  pm=pm + theme(axis.text.x = element_text(angle = 90, hjust = 1))

  return(pm)

}

### Flow Duration Curves (Appendix A4 b)

plotFDC=function(subsample, line_width=2, directory="~"){

  df=read.csv(here('case study data','hydrology_all_annual.csv'),
              stringsAsFactors = F)
  df$ID=paste(df$Scenario, df$TraceNumber, sep='.') # same as above, just for the original data frame

  all.filter=unique(df$ID)
  
  scatter.df=list()
  id.filter=list()
  q.ann.sample=list()
  for (i in 1:length(subsample)){
    scatter.df[[i]]=subsample[[i]]$SOW
    scatter.df[[i]]$ID=paste(scatter.df[[i]]$Scenario, scatter.df[[i]]$TraceNumber, sep='.')
    
    id.filter[[i]]=unique(scatter.df[[i]]$ID)
    q.ann.sample[[i]]=dplyr::filter(df, ID %in% id.filter[[i]]) # obtain the annual streamflow values from df contained in your efficient sample
    
    
  }
  
  all.filter=unique(df$ID)
  
  trans.const=.65
  full.col='gray'
  samp.col='deepskyblue'
  lwd=line_width
  
  
  c=1
  for (i in all.filter){ # loop through each unique scenario and trace
    filter.trace=dplyr::filter(df,ID==i)
    eqn=ecdf(filter.trace$CNF_LF)
    
    if(c==1){ #establish plot on first loop
      baseplot=plot(sort.default(filter.trace$CNF_LF), eqn(sort.default(filter.trace$CNF_LF)), xlab='Flow [MAF]', ylab = 'CDF',
                    main = 'Annual Flow at Lees Ferry', type = 'l', col=transparent(orig.col = full.col, trans.val =(trans.const+.1)), lwd=lwd)
    }
    else{ # add trace on subsequent loops
      lines(sort.default(filter.trace$CNF_LF), eqn(sort.default(filter.trace$CNF_LF)), xlab='Flow [MAF]', ylab = 'CDF',
            main = 'Annual Flow at Lees Ferry', type = 'l', col=transparent(orig.col = full.col, trans.val = (trans.const+.1)), lwd=lwd)
      
    }  
    c=c+1
  }
  
  baseplot=recordPlot() # record the full CDF as a baseplot
  
  #FDCplots=list()
  
  for (s in 1:length(subsample)){
    
    baseplot
    
    for(i in id.filter[[s]]){ # loop through efficient sample, add to plot
      filter.trace=dplyr::filter(df,ID==i)
      eqn=ecdf(filter.trace$CNF_LF)
      
      lines(sort.default(filter.trace$CNF_LF), eqn(sort.default(filter.trace$CNF_LF)), xlab='Flow [MAF]', ylab = 'CDF',
            main = 'Annual Flow at Lees Ferry', type = 'l', col=transparent(orig.col = samp.col, trans.val = trans.const), lwd=lwd)
      
    }
    
    legend('bottomright', legend = c('All traces', 'Sampled subset'),
           col = c(full.col,samp.col), 
           lty = c('solid', 'solid'), lwd=c(lwd,lwd),
           cex=0.89)
    
    #savePlot(filename = paste0(directory,"/plot",s, ".pdf"), type = "pdf")
    #savePlot(filename = paste0(directory,"/plot",s, ".jpeg"), type = "jpeg", device = )
    dev.off()
    
    #FDCplots[[s]]=recordPlot()
    
    
  }

  #return(FDCplots)
  
}


#### wrapper to plot_ly parallel coordinates (Appendix A3, A5)

par_coords=function(data, max_cols=NULL, n_var, color_var, title='User selected metrics', labels=colnames(data), source=NULL, policy_ID=NULL,
                    color_scale="blue-red", reverse_scale=F, legendTF=F, colorbarTitle="", maintainAxesRange=T, axes_data=data, show_ID_front=F, 
                    labelangle=-10, titlesize=20, labelsize=11, img_height=400, img_width=600, img_scale=10,
                    img_format="jpeg", categorical_cols=NULL, categorical_names_list=NULL){
  
  # identify maximization axes
  
  T_F=rep(FALSE,ncol(data))
  
  if (is.null(max_cols)){
    
  } else {
    
    T_F[max_cols]=TRUE
    
  }
  
  # identify categorical columns
  
  categorical=rep(FALSE, ncol(data))
  
  if(is.null(categorical_cols)){} else {categorical[categorical_cols]=T}
  

  # custom color palette
  mypal=matrix(nrow=4, ncol=2)
  mypal[,1]=c(0,.25,.75,1)
  mypal[,2]=c('rgb(181,177,177)','rgb(255,102,255)','rgb(153,255,255)','rgb(255,0,0)')
  
  
  dimensions=list()
  c=1 # counter to index categorical_names_list
  for (i in 1: ncol(data[,1:n_var])){
    
    dimensions[[i]]=list()
    
    if (maintainAxesRange==FALSE){ # use axes range of data. This means if policies are removed, axes ranges can shrink

      dimensions[[i]][['range']]=sort(range(data[[i]]), decreasing=T_F[i])
      
      
      if(categorical[i]){
        dimensions[[i]][["ticktext"]]=categorical_names_list[[c]]
        dimensions[[i]][["tickvals"]]=sort(unique(data[[i]]))
        c=c+1
      }
      
      
    } else { # axes ranges are kept as the range of all policies
     
      dimensions[[i]][['range']]=sort(range(axes_data[[i]]), decreasing=T_F[i])
      
      if(categorical[i]){
        dimensions[[i]][["ticktext"]]=categorical_names_list[[c]]
        dimensions[[i]][["tickvals"]]=sort(unique(axes_data[[i]]))
        c=c+1
      }
      
    }
    
    dimensions[[i]][['label']]=labels[i]
    dimensions[[i]][['values']]=data[[i]]
    
    if (colnames(data)[i]=='ID'){
      dimensions[[i]][['constraintrange']]=policy_ID
    }
    
    if (colnames(data)[i] %in% c("ID", "front")){
      dimensions[[i]][["visible"]]=show_ID_front
    }
    
    
  }
  
  
  if (color_scale=='blue-red'){color_scale=mypal}
  
  p <-data %>% plot_ly(type = 'parcoords', tickfont=list(size=13),
                       line = list(color =data[[color_var]], colorbar=list(title=list(text=colorbarTitle, side='right'), thickness=20, x=1.00, xpad=10),
                                   colorscale = color_scale, reversescale=reverse_scale, showscale=legendTF), labelangle=labelangle, labelside="top",
                       dimensions = dimensions, source=source, labelfont=list(size=labelsize)
                       
  ) # end plot_ly()
  
  # add title and set margins
  p=p %>% layout(margin=list(l=60,r=20,b=25,t=0, pad=0), title=list(text=title, font=list(size=titlesize)))
  
  # see these sources for ways to change the download image behavior. You can change height, width, and fyle type to obtain MUCH IMPROVED quality
  # https://plotly.com/python/configuration-options/
  # https://www.rdocumentation.org/packages/plotly/versions/4.9.3/topics/config
  # https://github.com/plotly/plotly.js/blob/master/src/plot_api/plot_config.js
  
  p=config(p, toImageButtonOptions=list(format=img_format, height=img_height, width=img_width, scale=img_scale))
  
  
  
  return(p)
  
  
}

