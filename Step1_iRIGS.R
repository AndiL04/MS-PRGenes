####Bayesian Framework for gene prediction
####Multiple Sclerosis
####Origin: Wang, Q., Chen, R., Cheng, F., Wei, Q., Ji, Y., Yang, H., Zhong, X., Tao, R., Wen, Z., Sutcliffe, J. S., Liu, C., Cook, E. H., Cox, N. J., & Li, B. (2019, May). 
####A Bayesian framework that integrates multi-omics data and gene networks predicts risk genes from schizophrenia GWAS data. Nat Neurosci, 22(5), 691-699. 
####https://doi.org/10.1038/s41593-019-0382-7

####Modified by: Andi Liu
####Last update date: 2/8/2021

###-----------------------------First Step-----------------------------------###
###Extract candidate genes based on GWAS information of MS###

####------------------------------GWAS to GENES----------------------------#####
###locate the address
setwd("C:/Users/andya/OneDrive/Desktop/Work/MS_Bayesian_submitted/iRIGS");

###read data
library(readxl)
gwas <- read_excel("aav7188_Patsopoulos_Tables S1-S10_S8.xlsx", 
                   sheet = "ST8", range = "A4:J204")  ###read 200 loci

colnames(gwas)[1:3] <- c("SNP","Chr","Pos_hg19")
gwas <- gwas[,c("SNP","Chr","Pos_hg19")]

###adding chr infront of chr number.
if(substr(gwas$Chr[1],1,3)!="chr")  gwas$Chr<-paste("chr",gwas$Chr,sep="")


###read gene infomation
gene <- read.delim("./supporting_files/All_human_genes")


###setting flanking time
flanking <- 1000000  ### Will find candidate genes within 1000000 basepair range

###extract candidate genes
output<-list()
for(i in 1:nrow(gwas))
{
  start<-max(gwas[i,]$Pos_hg19-flanking,0)
  end<-gwas[i,]$Pos_hg19+flanking
  temp<-gene[gene$chrom==gwas[i,]$Chr,]
  index1<-temp$start_hg19<end & temp$start_hg19>start
  index2<-temp$end_hg19<end & temp$end_hg19>start
  temp<-temp[index1|index2,]  
  
  if(nrow(temp)>0)
  {
    temp$SNP<-gwas[i,]$SNP
    temp$SNP_chr<-gwas[i,]$Chr
    temp$SNP_pos_hg19<-gwas[i,]$Pos_hg19
    output<-rbind(output,temp)
  }
} 
gene<-output
gene$Name<-substr(gene$Name,1,15)

###Result of extract candidate genes
transp<-"./supporting_files/go_propogation_probality_rp_0.3.RData" 
cat("Loading propagation probabilities...\n")
load(transp); 
nodes<-colnames(pro_p)


###Working on extracted candidate genes
gene<-gene[(!is.na(gene$official_name)),]
gene<-gene[is.element(gene$official_name,nodes),]
region<-split(gene$official_name,gene$SNP)
cat(paste(length(unique(gene$official_name)),"genes from",length(region),"loci were found with propagation probability...\n"))

pro_p<-pro_p[,(is.element(nodes,unique(gene$official_name)))]


####output candidate genes
bbb <- unique(gene$Name)
try <- paste(" -e ",bbb, sep = "")
write.table(t(try), "Unique candidate risk genes.txt", row.names = F, quote = F,col.names = F, sep = "")

####-----------------------END OF GWAS to GENES----------------------------#####
###-----------------------End of First Step----------------------------------###




###-----------------------------Second Step----------------------------------###
###------------------------Working with Covariates---------------------------###

####---------------------------- Distance to TSS ---------------------------####
compute_dist<-function(gene)
{
  index<-gene$strand=="+"
  tss_dist<-0
  tss_dist[index]<-abs(gene[index,]$start_hg19-gene[index,]$SNP_pos_hg19)
  tss_dist[!index]<-abs(gene[!index,]$end_hg19-gene[!index,]$SNP_pos_hg19)
  gene$dist<-tss_dist
  gene 
}

gene<-compute_dist(gene) 
####-----------------------END OF Distance to TSS --------------------------####

###---------------------------- adding FANTOM5  -----------------------------###

cat("Collecting extra evidence: FANTOM5...\n")

adding_fantom5<-function(gene)
{
  path<-"./supporting_files/Fantom5/" 
  eep<-paste(path,"enhancer_tss_associations.bed.txt",sep="")
  eep<-read.delim(eep,as.is=T)
  
  eep<-eep[unlist(lapply(strsplit(eep$name,";"),function(x) length(x)))==5,]
  pos<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[1]]))
  chr<-unlist(lapply(strsplit(pos,":"),function(x)x[[1]]))
  start<-as.numeric(lapply(strsplit(unlist(lapply(strsplit(pos,":"),function(x)x[[2]])),"-"),function(x) x[[1]]))
  end<-as.numeric(lapply(strsplit(unlist(lapply(strsplit(pos,":"),function(x)x[[2]])),"-"),function(x) x[[2]]))
  fdr<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[length(x)]]))
  fdr<-as.numeric(lapply(strsplit(fdr,":"),function(x) x[[2]]))
  genes<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[length(x)-2]]))
  
  eep$chr<-chr; eep$start<-start; eep$end<-end; eep$fdr<-fdr; eep$gene<-genes
  eep<-eep[,13:17];  eep<-eep[eep$fdr<1,]
  eep$enhancer<-paste(eep$chr,eep$start,eep$end,sep=":")
  
  gene$fantom5_enhancer_no<-0 
  for(i in 1:nrow(gene))
  {
    temp<-eep[eep$gene==gene[i,]$official_name,]
    
    if(nrow(temp)>0)
      gene[i,]$fantom5_enhancer_no<-length(unique(temp$enhancer))
  }
  gene
} 

gene<-adding_fantom5(gene)

###------------------------- end of adding FANTOM5  -------------------------###

### -------------------------   adding capture Hi-C  -----------------------####

cat("Collecting extra evidence: capture Hi-C...\n")

adding_caphic<-function(gene)
{
  path<-"./supporting_files/capHiC/"
  cap4<-paste(path,"GM12878_DRE_number",sep="")
  cap4<-read.delim(cap4,as.is=T)
  
  gene$cap4_enhancer_no<-cap4[match(gene$Name,cap4$Name),]$cap4_enhancer_no
  gene
}

gene<-adding_caphic(gene)

###----------------------- end of adding capture Hi-C  ----------------------###

###--------------------------- adding brain HiC -----------------------------###

cat("Collecting extra evidence: Brain Hi-C...\n")

adding_brainhic<-function(gene)
{
  path1<-"./supporting_files/BrainHiC/S22_TSS_CP.txt"
  path2<-"./supporting_files/BrainHiC/S23_TSS_GZ.txt"
  
  cp<-read.delim(path1,as.is=T)
  gz<-read.delim(path2,as.is=T)
  cp$enhancer<-paste(cp$chr,cp$interacting_bin_start,cp$interacting_bin_end,sep=":")
  gz$enhancer<-paste(gz$chr,gz$interacting_bin_start,gz$interacting_bin_end,sep=":")
  
  gene$brain_cp<-0; gene$brain_gz<-0 
  
  for(i in 1:nrow(gene))
  {
    temp<-cp[cp$ENSGID_for_TSS==gene[i,]$Name,]
    if(nrow(temp)>0)
      gene[i,]$brain_cp<-length(unique(temp$enhancer))
  }
  
  for(i in 1:nrow(gene))
  {
    temp<-gz[gz$ENSGID_for_TSS==gene[i,]$Name,]
    if(nrow(temp)>0)
      gene[i,]$brain_gz<-length(unique(temp$enhancer))
  }
  gene
}

gene<-adding_brainhic(gene)

###------------------------- end of adding brain HiC ------------------------###

###---------------------- adding differential expression --------------------###

dexpr<-read.delim("./supporting_files/MS_DE.csv",
                  sep = ",",row.names = 1)

###------------------ end of adding differential expression -----------------###


###---------------------- adding differential methylation --------------------###
methy <- read.csv("./supporting_files/methylation.ms.csv",
                  row.names = 1)

###------------------ end of adding differential methylation -----------------###



###-------------------------   combining p-values   -------------------------###

#extra<-gene[,c("dist","cap4_enhancer_no","fantom5_enhancer_no","tis_expr","brain_cp","brain_gz")]
extra<-gene[,c("dist","cap4_enhancer_no","fantom5_enhancer_no","brain_cp","brain_gz")]
#extra<-gene[,c("cap4_enhancer_no","fantom5_enhancer_no")]

for(i in 1:ncol(extra))
  extra[,i][is.na(extra[,i])]<-median(extra[,i][!is.na(extra[,i])]) 

sigma<-cov(extra)
s<-svd(sigma)
s1<- s$u %*% diag((s$d)^0.5) %*% t(s$v)
s2<- s$u %*% diag((s$d)^(-0.5)) %*% t(s$v)

mu<-apply(extra,2,mean)
extra_tran<-t(s2%*%(apply(extra,1,function(x) x-mu)))

lower_tail<-c(T,F,F,F,F)
for(i in 1:ncol(extra_tran))
  if(!lower_tail[i]) 
    extra_tran[,i]<-(-extra_tran[,i])

extra_p<-apply(extra_tran,2,pnorm)

### adding differential expression p-value
extra_p<-cbind(extra_p,dexpr[match(gene$official_name,rownames(dexpr)),]$pvalue)

### adding methylation p-value
extra_p<-cbind(extra_p,methy[match(gene$official_name,methy$gene),]$p_value)

for(i in 1:ncol(extra_p))
  extra_p[,i][is.na(extra_p[,i])]<-median(extra_p[,i][!is.na(extra_p[,i])])

gene$extra_weight<- (-log(apply(extra_p,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x),lower.tail=F)})))


extra_weight<-gene[,c("official_name","extra_weight")]
extra_weight<-extra_weight[match( unique(extra_weight$official_name),extra_weight$official_name),]
colnames(extra_weight)[1]<-"gene"
extra_weight
####------------------------- end of loading data ------------------------- ####



####================= burn in step ======================####

#####Bayesian Sampling
burnin_round <- 6000          ### maximum sampling round allowed at burn in step PRESET 3000
after_burnin_round <- 6000    ### maximum sampling round allowed at post burn in step PRESET 3000
exclude_samegene <- T

####================= burn in step ======================####

#t0<-proc.time() 

thres<-0.01; pickup<-0; 
num_region<-length(region); circle<-1; chosen<-NULL

#sampling initial gene set
remaining<-unlist(lapply(region,function(x) sample(x,1)))

num0<-rep(0,sum(unlist(lapply(region,length))))

dif<-thres+1; dif_record<-NULL 
#while(dif>thres && circle<100)
while(dif>thres && circle<(burnin_round+1))
{
  pickup<-pickup%%num_region+1  
  if(pickup==1)
    if(!is.null(chosen))
    {
      num1<-NULL
      for(j in 1:length(region))
        num1<-c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x)))) 
      num1<-num1+num0
      if(circle>1)
      {
        freq0<-num0/(num_region*(circle-1))
        freq1<-num1/(num_region*circle)
        dif<-(sum((freq0-freq1)^2))^0.5
        if( circle%%50==0 )
        {
          cat("Burnin sampling, sampling circle:",circle,"\n")
        }
        #dif_record<-c(dif_record,dif)
      }
      num0<-num1; chosen<-NULL; circle<-circle+1
    }
  
  pickup_p<-pro_p[,is.element(colnames(pro_p),remaining[-pickup])]
  pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]  
  
  if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
  {                                                                            ### conditional genes, exclude the same genes 
    pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
    if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
  }
  
  if(is.null(dim(pickup_p))) { pickup_p<- 1; names(pickup_p)<-region[[pickup]] } else    ### when there is only one candiate gene  
  { pickup_p<-apply(pickup_p,1,sum); pickup_p<-extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight*pickup_p } 
  
  if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i]<-1/length(pickup_p)     ### avoid probabilit=0 
  
  remaining[pickup]<-sample(names(pickup_p),1,replace=T,prob=pickup_p)
  chosen<-rbind(chosen,remaining)  
}

#proc.time()-t0 

###===================== end of burn in step ===============================###

###======================= post-burnin step ===================================###
#t0<-proc.time()

pickup<-0; num_region<-length(region); circle<-1; chosen<-NULL
num0<-rep(0,sum(unlist(lapply(region,length))))

joi_dis<-matrix(0,nrow=nrow(gene),ncol=nrow(gene))
temp<-NULL;
for(j in 1:length(region))
  temp<-c(temp,paste(names(region[j]),region[[j]],sep="_"))
colnames(joi_dis)<-temp; rownames(joi_dis)<-temp

thres<-0.01; dif<-thres+1
while(dif>thres && circle<(after_burnin_round+1) )
{
  pickup<-pickup%%num_region+1
  if(pickup==1)
    if(!is.null(chosen))  
    {
      ###================================= calculate frequency =========================###
      num1<-NULL
      for(j in 1:length(region))
        num1<-c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x))))
      num1<-num1+num0
      if(circle>1)
      {
        freq0<-num0/(num_region*(circle-1))
        freq1<-num1/(num_region*circle)
        dif<-(sum((freq0-freq1)^2))^0.5
        if( circle%%50==0 )
        {
          cat("Post-burnin sampling, sampling circle:",circle,"\n")
        }
        #dif_record<-c(dif_record,dif)
      }
      num0<-num1; circle<-circle+1; chosen<-NULL
      ###============================= end of calculating frequency =======================###
    }
  
  pickup_p<-pro_p[,is.element(colnames(pro_p),remaining[-pickup])]
  pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]
  
  if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
  {                                                                            ### conditional genes, exclude the same genes
    pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
    if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
  }
  
  if(is.null(dim(pickup_p))) { pickup_p<-1; names(pickup_p)<-region[[pickup]] } else    ### when there is only one candiate gene
  { pickup_p<-apply(pickup_p,1,sum); pickup_p<-extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight*pickup_p } 
  
  if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i]<-1/length(pickup_p)     ### avoid probabilit=0
  
  remaining[pickup]<-sample(names(pickup_p),1,replace=T,prob=pickup_p)
  chosen<-rbind(chosen,remaining)
  
  ###=============================== calculating joint distribution  =====================###
  index_col<-match(paste(names(remaining[-pickup]),remaining[-pickup],sep="_"),colnames(joi_dis))  
  index_row<-match(paste(names(remaining[pickup]),remaining[pickup],sep="_"),colnames(joi_dis))  
  joi_dis[index_row,index_col]<-joi_dis[index_row,index_col]+1 
}

#proc.time()-t0

###=====================  end of post-burnin step  =======================================###

###===================== summarize and record the results ================================###
freq<-cbind(unlist(region),freq1)

region_indicator<-NULL
gene_num<-as.numeric(lapply(region,length))
for(i in 1:length(gene_num))
  region_indicator<-c(region_indicator,rep(names(region[i]),gene_num[i]))

freq<-cbind(freq,region_indicator)
colnames(freq)<-c("gene","post_prob","region")
freq<-as.data.frame(freq,stringsAsFactors=F)
freq[,2]<-as.numeric(freq[,2])

output<-NULL                                      #### sort according to posterior probability
for(i in unique(freq$region))
{
  temp<-freq[freq$region==i,]
  output<-rbind(output,temp[order(temp$post_prob,decreasing=T),])
}
freq<-output

#####output results####
cat("Recording the results!\n")
#freq <- read.delim("MS_GWAS_Risk_Genes_methy.txt")
write.table(freq,"MS_GWAS_Risk_Genes_methy.txt",quote=F,row.names=F,sep="\t")


###output genelist with highest posterior probabilities###
geneset <- data.frame()
for (i in unique(freq$region)){
  reg <- paste(i)
  m <- freq[freq$region == i,]
  geneset <- rbind(geneset,m[1,])
}
genelist <- unique(geneset$gene)

write.table(genelist, "MS_Risk_GeneList_methy.txt", quote = F, row.names = F, sep = "\t")
write.csv(geneset, "MS_risk_geneset.csv",row.names = F)


