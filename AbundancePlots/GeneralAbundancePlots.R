cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")

obj<-ps_filt

obj<-transform_sample_counts(obj, function(x) x / sum(x))

i<-which(colnames(tax_table(obj))==tax)
taxa_names(obj)<-tax_table(obj)[,colnames(tax_table(obj))==tax]

titl<-paste(tax,"Relative Abundance")
a<-plot_bar(obj,fill=tax,title = titl)
a<-a+theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())
print(a)

taxa_names(obj)<-tax_table(obj)[,colnames(tax_table(obj))==tax]

prop_data<-data.frame(otu_table(obj))
row.names(prop_data)<-tax_table(obj)[,colnames(tax_table(obj))==tax]

tab<-data.frame(t(prop_data),Experiment=sample_data(obj)$Experiment,
                IDs=rownames(sample_data(obj))) 

tab<-melt(tab, id=c("Experiment","IDs"))

b<-ggplot(data=tab,aes(y=variable,x=value))+
  geom_point(aes(color=Experiment),size=2)+
  xlab("Relative Abundance")+
  ylab("Taxa")+
  ggtitle(tax)+
  scale_colour_manual(values=cbPalette)+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom",legend.title = element_blank())

print(b)
