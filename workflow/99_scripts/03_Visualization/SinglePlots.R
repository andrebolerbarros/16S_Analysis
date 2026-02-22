criteria<-paste0("detected_",cut.off,sep="")
coll<-paste0("W_frame$",criteria,sep="")
j<-which(colnames(W_frame)==criteria)

tab<-data.frame(prop_data,Profile=sample_data(obj)$Profile,Treat_Surv=sample_data(obj)$Treat_Surv,
                IDs=rownames(sample_data(obj))) 

tab<-melt(tab, id=c("Profile","IDs","Treat_Surv"))

all(W_frame$otu_names %in% tab$variable)

temp<-W_frame[W_frame[(which(colnames(W_frame)==criteria))]==TRUE,]
temp<-tab[tab$variable %in% temp$otu_names,]

g1<-ggplot(data=temp,aes(y=value,x=Profile))+
  geom_point(aes(color=Treat_Surv),size=2)+
  xlab("Profile")+
  ylab("Relative Abundance")+
  scale_colour_manual(values=cbPalette)+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom",legend.title = element_blank())+
  facet_wrap( ~ variable,scales = "free")

print(g1)
