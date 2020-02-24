######################################################################################################################
R-script for plotting bootstrapped psmc trajectories for different t values 

Written By: Ajinkya Bharatraj Patil.
######################################################################################################################

library("scales")

t500 <- read.table("../t500/RES1386-IW.rg.qf.s100.t500.bp.g22.0.txt")
t65 <- read.table("../t65/RES1386-IW.rg.qf.s100.t65.g22.0.txt")
t55 <- read.table("../t55/RES1386-IW.rg.qf.s100.t55.g22.0.txt")
t45 <- read.table("../t45/RES1386-IW.rg.qf.s100.t45.g22.0.txt")
t35 <- read.table("../t35/RES1386-IW.rg.qf.s100.t35.g22.0.txt")

t500_bs = list.files("../t500",pattern="RES1386-IW.rg.qf.s100.t500.bp.combined.g22.*.txt", full.names = TRUE)
t65_bs = list.files("../t65",pattern="RES1386-IW.rg.qf.s100.t65.combined.g22.*.txt", full.names = TRUE)
t55_bs = list.files("../t55",pattern="RES1386-IW.rg.qf.s100.t55.combined.g22.*.txt", full.names = TRUE)
t45_bs = list.files("../t45",pattern="RES1386-IW.rg.qf.s100.t45.combined.g22.*.txt", full.names = TRUE)
t35_bs = list.files("../t35",pattern="RES1386-IW.rg.qf.s100.t35.combined.g22.*.txt", full.names = TRUE)

pdf("COALT2_BS.pdf",height=10,width=15)
par(mar = c(6, 4, 4, 4))
plot(t500$V1/1000,t500$V2,type="n",log="x",lwd=3,xaxt="n",yaxt="n",ylim=c(0,10),main="Effect of variable TMRCA values in PSMC on Mesua fe
rrea genome",xlab="",ylab="",axes=FALSE,cex.main=1.75)

rect(0.1,0,12.7,10, col=scales::alpha(rgb(col2rgb("cyan")[1,],col2rgb("cyan")[2,],col2rgb("cyan")[3,],max = 255), 0.25),border="transpare
nt")
rect(12.8,0,2580,10, col=scales::alpha(rgb(col2rgb("yellow")[1,],col2rgb("yellow")[2,],col2rgb("yellow")[3,],max = 255), 0.25),border="tr
ansparent")
rect(2581,0,7550,10, col=scales::alpha(rgb(col2rgb("green")[1,],col2rgb("green")[2,],col2rgb("green")[3,],max = 255), 0.25),border="trans
parent")
rect(11.7,0,115,10, col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],max = 255), 0.25),border="transparent"
)
rect(374,0,424,10, col = "gray", border = "transparent")
rect(300,0,337,10, col = "gray", border = "transparent")
rect(191,0,243,10, col = "gray", border = "transparent")
rect(82,0,123,10, col = "gray", border = "transparent")
text(100,9.2,"MIS 5",cex=1,srt=90,font=2)
text(4,9.2,"Holocene",cex=1,font=2)
text(30,9.2,"Last Glacial Period",cex=1,font=2)
text(2580/2,9.2,"Pleistocene",cex=1,font=2)
text(3300,9.2,"Pliocene",cex=1,srt=90,font=2)
text(400,9.2,"MIS 11",cex=1,srt=90,font=2)
text(320,9.2,"MIS 9",cex=1,srt=90,font=2)
text(225,9.2,"MIS 7",cex=1,srt=90,font=2)

for (i in t500_bs){
read.table(file=i,header=FALSE)-> t500_bs_all
lines(t500_bs_all$V1/1000,t500_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],
max = 255), 0.25))
}
lines(t500$V1/1000,t500$V2,type="s",col="red",lwd=3)

for (i in t65_bs){
read.table(file=i,header=FALSE)-> t65_bs_all
lines(t65_bs_all$V1/1000,t65_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("brown")[1,],col2rgb("brown")[2,],col2rgb("brown")[3,],max = 255), 0.25))
}
lines(t65$V1/1000,t65$V2,type="s",col="brown",lwd=3)

for (i in t55_bs){
read.table(file=i,header=FALSE)-> t55_bs_all
lines(t55_bs_all$V1/1000,t55_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("green")[1,],col2rgb("green")[2,],col2rgb("green")[3,],max = 255), 0.25))
}
lines(t55$V1/1000,t55$V2,type="s",col="green",lwd=3)

for (i in t45_bs){
read.table(file=i,header=FALSE)-> t45_bs_all
lines(t45_bs_all$V1/1000,t45_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("blue")[1,],col2rgb("blue")[2,],col2rgb("blue")[3,],max = 255), 0.25))
}
lines(t45$V1/1000,t45$V2,type="s",col="blue",lwd=3)

for (i in t35_bs){
read.table(file=i,header=FALSE)-> t35_bs_all
lines(t35_bs_all$V1/1000,t35_bs_all$V2, lty=3, type="s",col=scales::alpha(rgb(col2rgb("cyan")[1,],col2rgb("cyan")[2,],col2rgb("cyan")[3,],max = 255), 0.25))
}
lines(t35$V1/1000,t35$V2,type="s",col="cyan",lwd=3)

legend("topleft",legend=c("-t 500","-t 65", "-t 55", "-t 45", "-t 35"),fill=c("red","brown","green","blue","cyan"),cex=1.5)

axis(1, col="black",las=1,font=2,font.axis=2)
mtext(side=1, line=2, "Time in past (KYA)", col="black", font=2,cex=1.25)
axis(2, col="black",las=1,font=2,font.axis=2)
mtext(side=2, line=2, "Effective population size (Ne) x 10e4",col="black", font=2,cex=1.25)
box()
dev.off()
