library("scales")

hg38 <- read.table("../bgi_illumina_bs/hg38.0.txt")
bgi <- read.table("../bgi_illumina_bs/bgi_human.un.0.txt")

bgi_bs = list.files("../bgi_illumina_bs/",pattern="bgi_human.combined.*.txt", full.names = TRUE)
hg38_bs = list.files("../bgi_illumina_bs/",pattern="hg38.combined.*.txt", full.names = TRUE)


pdf("BGI_illumina_BS.pdf",height=10,width=20)
par(mar = c(6, 4, 9, 6))
plot(hg38$V1,hg38$V2,type="n",log="x",lwd=3,ylim=c(0,4), xaxt="n",yaxt="n",main="BGISEQ vs ILLUMINA effect on Human Genome Demographic Inference",xlab="",ylab="",axes=FALSE,cex.main=1.75)

for (i in bgi_bs){
read.table(file=i,header=FALSE)-> bgi_bs_all
lines(bgi_bs_all$V1,bgi_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],max = 255), 0.25))
}
lines(bgi$V1,bgi$V2,type="s",col="red",lwd=3)

for (i in hg38_bs){
read.table(file=i,header=FALSE)-> hg38_bs_all
lines(hg38_bs_all$V1,hg38_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("blue")[1,],col2rgb("blue")[2,],col2rgb("blue")[3,],max = 255), 0.25))
}
lines(hg38$V1,hg38$V2,type="s",col="blue",lwd=3)

legend("top",legend=c("hg38-Illumina","hg38-BGISEQ"),fill=c("blue","red"),cex=1.5)

axis(1, col="black",las=1,font=2,font.axis=2)
mtext(side=1, line=2, "Years ago", col="black", font=2,cex=1.25)
axis(2, col="black",las=1,font=2,font.axis=2)
mtext(side=2, line=2, "Effective population size (Ne) x 10e4",col="black", font=2,cex=1.25)
box()

dev.off()
