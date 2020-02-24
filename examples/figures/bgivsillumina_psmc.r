#script for plotting BGISEQ vs ILLUMINA PSMC plots
#Written By : Ajinkya Bharatraj Patil

bgi_t5 <- read.table("../bgi_human.t5.0.txt")
bgi_t15 <- read.table("../bgi_human.t15.0.txt")
ilm_t15 <- read.table("../SRR9091899.t15.0.txt")
ilm_t5 <- read.table("../SRR9091899.t5.0.txt")
pdf("BGI_v_ILLUMINA_PSMC.pdf",height=10,width=15)
plot(ilm_t15$V1,ilm_t15$V2,log="x",type="n",xlab="log(Time in past)",ylab="Effective Population size x 10e04",main="BGISEQ vs ILLUMINA PSMC")
lines(ilm_t15$V1,ilm_t15$V2,type="s",col="red")
lines(ilm_t5$V1,ilm_t5$V2,type="s",col="orangered")
lines(bgi_t15$V1,bgi_t15$V2,type="s",col="blue")
lines(bgi_t5$V1,bgi_t5$V2,type="s",col="skyblue")
legend("topleft",legend=c("ILLUMINA_t5","ILLUMINA_t15","BGISEQ_t5","BGISEQ_t15"),fill=c("orangered","red","skyblue","blue"),cex=1.5)
dev.off()
