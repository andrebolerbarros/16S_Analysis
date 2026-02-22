library(ggplot2)
theme_set(theme_bw())
#######################
filt<-data.frame(Criteria = 0,Remaining=dim(otu_table(obj))[1])
colnames(filt)<-c("Criteria","Remaining")

for (i in seq(0:20)){
  j<-i+1
  a<-otu_table(obj)[which(taxa_sums(obj) / sum(taxa_sums(obj)) > i *1E-5), ]
  c<-dim(a)[1]
  filt[j,]<-c((i*1E-5)*100,c)
}

g1<-qplot(x = filt$Criteria,y=filt$Remaining,
          xlab = "Relative Abundance Threshold",ylab="Nr. Remaining OTU's",main = "Filtering Results based on Relative Abundance of OTUs (%)")+
  scale_y_continuous(breaks = seq(0,1400,200))+
  geom_hline(yintercept = c(500,400,300),linetype="dashed",color=c("limegreen","blue","red"))+
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5))+
  xlim(-0.001,0.021)

print(g1)

keepTaxa = otu_table(obj)[which(taxa_sums(obj) / sum(taxa_sums(obj)) > 5 *1E-5), ]

obj_filt<-prune_taxa(taxa_names(keepTaxa),obj)

print(paste0("Remaining OTU's: ", dim(otu_table(obj_filt))[1],sep=""))
