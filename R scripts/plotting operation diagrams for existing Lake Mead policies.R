# plotting the current Lake Mead operations
# Nathan Bonham
# March 2022


# Combined operation = Interim Guidelines + BWSCP + DCP

policy_lab=c("07 IG" ,"Min 323","DCP", "BWSCP", "Combined")
policy_lab_y=c(1075,1075,1090,1090,1090)
T1e=policy_lab_y
maxVol=c(500, 125, 600, 150, 1375)
T1V=c(333, 50, 200, 41, 241)
nTiers=c(3,3,5,8,8)
ID=1:5

existing_wide=data.frame(policy_lab, policy_lab_y, T1e,maxVol, T1V,nTiers,ID)

policy=as.character(rep(c(1:5), each=9))
Tier=rep(c(1,2,3,4,5,6,7,8, "dead pool"), 5)

IGe=c(1075,1050,1025,895,895,895,895,895,895)
Min323e=c(1075,1050,1025,895,895,895,895,895,895)
DCPe=c(1090,1045,1040,1035,1030,895,895,895,895)
BWSCPe=c(1090,1075,1050,1045,1040,1035,1030,1025, 895)
Combinede=c(1090,1075,1050,1045,1040,1035,1030,1025, 895)


IDdelta=IGe[-9]-IGe[-1]
Min323delta=Min323e[-9]-Min323e[-1]
DCPdelta=DCPe[-9]-DCPe[-1]
BWSCPdelta=BWSCPe[-9]-BWSCPe[-1]
Combineddelta=Combinede[-9]-Combinede[-1]

delta=c(IDdelta,895, Min323delta, 895,DCPdelta,895, BWSCPdelta, 895,Combineddelta,895)

IGv=c(333,417,500,0,0,0,0,0,0)
Min323v=c(50, 70, 125, 0,0,0,0,0,0)
DCPv=c(200,450, 500, 550, 600, 0,0,0,0) # this includes Minute 323, BWSCP, and DCP
BWSCPv=c(41,30,34,76,84,92,101,150,0)
Combinedv=c(241,613,721,1013,1071,1129,1188,1375,0)

elevation=c(
  IGe,
  Min323e,
  DCPe,
  BWSCPe,
  Combinede
)

volume=c(IGv, Min323v,DCPv, BWSCPv,Combinedv)

existing_long=data.frame(policy,Tier, delta, v_lab=paste (volume, ""), elevation, volume)

existing_policies=DV_plot(long.data = existing_long, wide.data = existing_wide,to_plot = 1:5,metric = NULL, xlabel = "",ID_label = NULL,y_axis2 = F,shrink_legend = F, interactive = F, volume_labs = T, labelsize = 4.5, hex_shift = F)
existing_policies=existing_policies+ylab("pool elevation (ft msl)")
existing_policies

# plot the combined policy compared to 3 very similar policies in neuron 3

similar=c(85, 233, 416) # IDs of policies in neuron 3 with very similar T1e, T1V, maxV compared to combined

compare_long=rbind(dplyr::filter(existing_long, policy==5), dplyr::filter(long_data, policy %in% similar))
temp=str_split(string = compare_long$v_lab, pattern = " ") # creates list
compare_long$v_lab=do.call(rbind.data.frame, temp)[,1] # convert to df, add to compare_long
compare_wide=rbind(dplyr::filter(existing_wide, ID==5), dplyr::filter(wide_data, ID %in% similar)[,-2])
compare_policies=DV_plot(long.data = compare_long, wide.data = compare_wide,to_plot = c(5,similar),metric = NULL, xlabel = "",
                         ID_label = NULL,y_axis2 = F,shrink_legend = F, interactive = F,
                         volume_labs = T,v_lab_nudge = -2, labelsize = 3, hex_shift = F, summary_stats = F )
compare_policies=compare_policies+ylab("pool elevation (ft msl)")
compare_policies



