#script for plotting bootstrapped human genome assembly based PSMC
#Written By: Ajinkya Bharatraj Patil.

library("scales")
hg4 <- read.table("../hg4/hg4.0.txt")
hg10 <- read.table("../hg10/hg10.0.txt")
hg15 <- read.table("../hg15/hg15.0.txt")
hg19 <- read.table("../hg19/hg19.0.txt")
hg38 <- read.table("../bgi_illumina_bs/hg38.0.txt")

hg4_bs = list.files("../hg4",pattern="hg4.combined.*.txt", full.names = TRUE)
hg10_bs = list.files("../hg10",pattern="hg10.combined.*.txt", full.names = TRUE)
hg15_bs = list.files("../hg15",pattern="hg15.combined.*.txt", full.names = TRUE)
hg19_bs = list.files("../hg19",pattern="hg19.combined.*.txt", full.names = TRUE)
hg38_bs = list.files("../bgi_illumina_bs",pattern="hg38.combined.*.txt", full.names = TRUE)

pdf("Human_assmblies_BS.pdf",height=10,width=20)
par(mar = c(6, 4, 9, 6))
plot(hg4$V1,hg4$V2,type="n",log="x",lwd=3,ylim=c(0,4), xaxt="n",yaxt="n",main="Effect of Human Genome Assembly Quality on Demographic Inf
erence",xlab="",ylab="",axes=FALSE,cex.main=1.75)

for (i in hg4_bs){
read.table(file=i,header=FALSE)-> hg4_bs_all
lines(hg4_bs_all$V1,hg4_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],max = 255), 0.25))
}
lines(hg4$V1,hg4$V2,type="s",col="red",lwd=3)

for (i in hg10_bs){
read.table(file=i,header=FALSE)-> hg10_bs_all
lines(hg10_bs_all$V1,hg10_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],max = 255), 0.25))
}
lines(hg10$V1,hg10$V2,type="s",col="red",lwd=3)

for (i in hg15_bs){
read.table(file=i,header=FALSE)-> hg15_bs_all
lines(hg15_bs_all$V1,hg15_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],max = 255), 0.25))
}
lines(hg15$V1,hg15$V2,type="s",col="red",lwd=3)

for (i in hg19_bs){
read.table(file=i,header=FALSE)-> hg19_bs_all
lines(hg19_bs_all$V1,hg19_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],max = 255), 0.25))
}
lines(hg19$V1,hg19$V2,type="s",col="red",lwd=3)

for (i in hg38_bs){
read.table(file=i,header=FALSE)-> hg38_bs_all
lines(hg38_bs_all$V1,hg38_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("blue")[1,],col2rgb("blue")[2,],col2rgb("blue")[3,],max = 255), 0.25))
}
lines(hg38$V1,hg38$V2,type="s",col="blue",lwd=3)

legend("top",legend=c("hg4 (92.86%)","hg10 (97.44%)","hg15 (99.84%)","hg19 (100%)","hg38 (100%)"),fill=c("red","orangered","plum1","purpl
e","blue"),cex=1.5)

axis(1, col="black",las=1,font=2,font.axis=2)
mtext(side=1, line=2, "Years ago", col="black", font=2,cex=1.25)
axis(2, col="black",las=1,font=2,font.axis=2)
mtext(side=2, line=2, "Effective population size (Ne) x 10e4",col="black", font=2,cex=1.25)
box()
par(new=TRUE)
plot(c(1:58),c(1:58),ylim=c(0,0.3),type="n",xlab="",ylab="",axes=FALSE,font=2,font.axis=2)
rep_col <- c("red","blue")
colcount <- 0
colcount <- colcount+1
read.table(file="../hg4/hg4.agg")->A
as.data.frame(aggregate(A$V2,list(A$V6),sd))->B
as.data.frame(aggregate(A$V2,list(A$V6),mean))->Bm
B$x/Bm$x->B$cv
lines(B$Group.1,B$cv,col=rep_col[colcount])
colcount <- colcount+1
read.table(file="../hg10/hg10.agg")->A
as.data.frame(aggregate(A$V2,list(A$V6),sd))->B
as.data.frame(aggregate(A$V2,list(A$V6),mean))->Bm
B$x/Bm$x->B$cv
lines(B$Group.1,B$cv,col=rep_col[colcount])
colcount <- colcount+1
read.table(file="../hg15/hg15.agg")->A
as.data.frame(aggregate(A$V2,list(A$V6),sd))->B
as.data.frame(aggregate(A$V2,list(A$V6),mean))->Bm
B$x/Bm$x->B$cv
lines(B$Group.1,B$cv,col=rep_col[colcount])
colcount <- colcount+1
read.table(file="../hg19/hg19.agg")->A
as.data.frame(aggregate(A$V2,list(A$V6),sd))->B
as.data.frame(aggregate(A$V2,list(A$V6),mean))->Bm
B$x/Bm$x->B$cv
lines(B$Group.1,B$cv,col=rep_col[colcount])
colcount <- colcount+1
read.table(file="../bgi_illumina_bs/hg38.agg")->A
as.data.frame(aggregate(A$V2,list(A$V6),sd))->B
as.data.frame(aggregate(A$V2,list(A$V6),mean))->Bm
B$x/Bm$x->B$cv
lines(B$Group.1,B$cv,col=rep_col[colcount])
mtext("Atomic Interval",side=3,line=2.5,font=2,font.axis=2,cex=1.25)
mtext("CV across Ne of bootstraps",side=4,line=2.5,font=2,font.axis=2,cex=1.25)
axis(3, ylim=c(1,58),las=1,font=2,font.axis=2)
axis(4, ylim=c(1,58),las=1,font=2,font.axis=2)

dev.off()
