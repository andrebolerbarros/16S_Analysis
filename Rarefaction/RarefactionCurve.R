samples<-vector()
obs<-vector()
sh<-vector()
depth<-vector()


for (i in seq(1,max(sample_sums(obj)),by = 500)) {
  rarified_obj <- rarefy_even_depth(obj, sample.size = i,replace = F,verbose = F,rngseed = 1)
  alpha_diversity <- estimate_richness(rarified_obj, measures = c("Observed","Shannon"))
  size<-length(alpha_diversity$Observed)
  
  alpha_diversity$depth<-rep(i,times=size)
  
  samples<-c(samples,rownames(alpha_diversity))
  obs<-c(obs,alpha_diversity$Observed)
  sh<-c(sh,alpha_diversity$Shannon)
  depth<-c(depth,rep(i,size))
}

tab_obs<-data.frame(depth,samples,obs)
tab_sh<-data.frame(depth,samples,sh)

g1<-ggplot(tab_obs,aes(x=depth,y=obs,colour=samples,group=samples))+
  geom_line(size = 1)+
  xlab("Nr. of Sequences")+
  ylab("Observed OTUs")+
  scale_x_continuous(breaks = seq(0,60000,10000))+
  geom_vline(xintercept = min(sample_sums(obj)),linetype="dashed",color="blue")+
  ggtitle("Rarefaction Curves (Observed OTUs)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "none")

g2<-ggplot(tab_sh,aes(x=depth,y=sh,colour=samples,group=samples))+
  geom_line(size = 1)+
  xlab("Nr. of Sequences")+
  ylab("Shannon Index")+
  scale_x_continuous(breaks = seq(0,60000,10000))+
  geom_vline(xintercept = min(sample_sums(obj)),linetype="dashed",color="blue")+
  ggtitle("Rarefaction Curves (Shannon Index)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "none")

grid.arrange(g1,g2,nrow=2)
