#根据差异基因在keggMap中进行着色
#差异上调为红色，差异下调为绿色
library(KEGGREST)
library(pathview)
source("C:/Users/asus/Documents/R/win-library/3.5/pathview/R/my_pathview.R")

oga = "ath"    #物种
path = 'C:\\Users\\asus\\Desktop\\1\\kegg'     #工作路径
setwd(path)

#data是差异基因表格，必须有AccID、Style、KeggID三列，即基因名、基因上下调信息和基因对应的keggID信息
data = read.table("5Avs5B_DEGs.ath.txt",sep="\t",header=T, fill=TRUE, quote=NULL)
cpt = read.table("cpt.txt",sep="\t",header=T, fill=TRUE, quote=NULL)  #代谢物表格，第一列是代谢物名字，第二列是数值
cpd=as.data.frame(cpt[,-1])  #代谢物格式要转成named num
cpd =t(cpd)
cpd = as.numeric(cpd)
names(cpd)=cpt[,1]

data2=data.frame(data$AccID,data$Style,data$KeggID)
colnames(data2)=c("AccID","Style","KeggID")

allPath=keggLink("pathway",oga)                #获取物种所有基因Pathway
allPath2=as.data.frame(allPath)
allPath2[,2]=names(allPath)
colnames(allPath2)=c("allpath","KeggID")

allko=keggLink("ko",oga)                #获取物种所有基因ko
allko2=as.data.frame(allko)
allko2[,2]=names(allko)
colnames(allko2)=c("ko","KeggID")

geneid=keggConv(oga,"ncbi-geneid")            #获取物种所有基因NCBI geneid
geneid2=as.data.frame(geneid)
geneid2[,2]=names(geneid)
idcov=function(x) {gsub('ncbi-geneid:','',x)}
geneid2=apply(geneid2,2,idcov)
colnames(geneid2)=c("KeggID","geneid")
geneid2=as.data.frame(geneid2)

path_geneid=merge(allPath2, geneid2, all=TRUE)
path_ko=merge(path_geneid, allko2, all=TRUE)
new_data=merge(data,path_ko,all.x=TRUE)
write.table(new_data,quote=FALSE, row.names=FALSE, sep="\t", eol="\n",,file="GeneInfo.txt")
data3=merge(data2,path_geneid,all.x=TRUE)
up <- dplyr::filter(data3, Style=="up")      #上调记为1，下调记为-1，即上调红色，下调绿色
up <- dplyr::mutate(up, style = 1)
down <- dplyr::filter(data3, Style=="down")
down <- dplyr::mutate(down, style = -1)
data3 <- rbind(up,down)

data4 <- dplyr::select(data3, geneid,style,allpath)
colnames(data4)=c("NCBIgeneid","Style","PathID")
data5=na.omit(data4)
pathID=as.data.frame(data5$PathID)
pathidcov=function(x) {gsub(paste("path:",oga,sep = ""),'',x)}
pathID2=apply(pathID,2,pathidcov)
pathID2=pathID2[!duplicated(pathID2)]

data6=data5[,-3]
data6=data6[!duplicated(data6),]
data7=as.data.frame(data6[,2])
row.names(data7)=data6[,1]
data7=as.matrix(data7)

for (i in 1:length(pathID2)){
  pv.out <- pathview(gene.data = data7[, 1], pathway.id =pathID2[i], cpd.data = cpd, species =oga, min.nnodes = 1,
                     out.suffix = "KeggMap", kegg.native =TRUE, plot.col.key=FALSE)
}

