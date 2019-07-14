library(exactRankTests)
library(nlme)
library(ggplot2)
###############################################
#####################################################
obj<-ps_filt

i<-which(colnames(tax_table(obj))==tax)
taxa_names(obj)<-tax_table(obj)[,colnames(tax_table(obj))==tax]

otu_data=data.frame(t(otu_table(obj))) #OTU list w/ sample ID
var_data=data.frame(sample_data(obj)) #metadata in fo

all(row.names(otu_data)==row.names(var_data)) #check if they match completely

n_otu=dim(otu_data)[2] #Number of OTU's, subtracting the "sample ID" col
otu_ids=colnames(otu_data) #OTU's ID, subtracting the "sample ID" col

logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu) #create a squared matrix for the logratios
for(ii in 1:(n_otu-1)){ #seq from 1 to Number of otu's - 1
  for(jj in (ii+1):n_otu){ #seq from 2 to Number of OTU's
    data.pair <- otu_data[,which(colnames(otu_data)%in%otu_ids[c(ii,jj)])] #create a data frame w/ the 2 OTU's
    data.pair$Profile<- var_data$Profile
    
    data.pair$lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2]))) #create a log-ratio column, with pseudo-counts of 1
    
    normality<-shapiro.test(data.pair$lr)
    #print(normality$p.value)
    
    if (normality$p.value > 0.05) {
      model=t.test(data.pair$lr ~ data.pair$Profile)
    } else {model = wilcox.test(data.pair$lr ~ data.pair$Profile)}
    
    logratio.mat[ii,jj] <- model$p.value #select the p-value from the row for the specific variable
  }}

ind <- lower.tri(logratio.mat)
logratio.mat[ind] <- t(logratio.mat)[ind]

logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1

mc.pval <- t(apply(logratio.mat,1,function(x){
  s <- p.adjust(x, method = "BH")
  return(s) }))

sig<-0.05

W <- apply(mc.pval,1,function(x){
  subp <- sum(x < sig,na.rm=T)
})

test.pval<-W > 0

p.zeros<-data.frame(1,1)
for (i in 1:ncol(otu_data)) {
  tab<-data.frame(taxa=otu_data[,i],info=var_data$Profile)
  tab<-aggregate(tab$taxa,by=list(tab$info), FUN=function(x){
    s=length(which(x==0))/length(x)})
  p.zeros[i,]<-tab$x
  rownames(p.zeros)[i]<-colnames(otu_data)[i]
}


av_counts<-data.frame(1,1)
eff_size<-vector()
prop<-transform_sample_counts(obj, function(x) x / sum(x))
prop_data<-t(otu_table(prop))
colnames(prop_data)<-tax_table(obj)[,colnames(tax_table(obj))==tax]

for (i in 1:ncol(otu_data)){
  tab<-data.frame(taxa=otu_data[,i],info=var_data$Profile)
  tab<-aggregate(tab[,1],by=list(tab$info),FUN=mean,na.rm=T)
  av_counts[i,]<-tab$x
  
  tab<-data.frame(taxa=prop_data[,i],info=var_data$Profile)
  tab<-aggregate(tab[,1],by=list(tab$info),FUN=median, na.rm=TRUE)
  eff_size[i]<-(tab$x[1]+1E-5)/(tab$x[2]+1E-5)
}


W_frame = data.frame(colnames(otu_data),av_counts,p.zeros,eff_size,W,row.names=NULL)
colnames(W_frame)<-c("otu_names","Average Counts Acla","Average Counts Epi","% Zeros Acla","% Zeros Epi","Effect Size","W_stat")

W_frame$detected_0.90=rep(FALSE,dim(W_frame)[1])
W_frame$detected_0.80=rep(FALSE,dim(W_frame)[1])
W_frame$detected_0.75=rep(FALSE,dim(W_frame)[1])
W_frame$detected_0.70=rep(FALSE,dim(W_frame)[1])

W_frame$detected_0.90[which(W_frame$W_stat>0.90*(length(colnames(otu_data))-1))]=TRUE
W_frame$detected_0.80[which(W_frame$W_stat>0.80*(length(colnames(otu_data))-1))]=TRUE
W_frame$detected_0.75[which(W_frame$W_stat>0.75*(length(colnames(otu_data))-1))]=TRUE
W_frame$detected_0.70[which(W_frame$W_stat>0.70*(length(colnames(otu_data))-1))]=TRUE

