#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
mu <- as.numeric(args[2])
gen <- as.numeric(args[3])
gen_size <- as.numeric(args[4])

#######################################################################################################################################
###################################### Reading in the input files #####################################################################
#######################################################################################################################################
pt_hm_s20 <- read.table(paste("s20/",args[1],".hm.s20.0.txt",sep=""))
pt_hm_s50 <- read.table(paste("s50/",args[1],".hm.s50.0.txt",sep=""))
pt_hm_s100 <- read.table(paste("s100/",args[1],".hm.s100.0.txt",sep=""))
pt_un_s20 <- read.table(paste("s20/",args[1],".un.s20.0.txt",sep=""))
pt_un_s50 <- read.table(paste("s50/",args[1],".un.s50.0.txt",sep=""))
pt_un_s100 <- read.table(paste("s100/",args[1],".un.s100.0.txt",sep=""))
msmc_hm <- read.table(paste("msmc/",args[1],".hm.final.txt",sep=""),header=T)
msmc_un <- read.table(paste("msmc/",args[1],".un.final.txt",sep=""),header=T)
hz_un_s100 <- read.table(paste("s100/",args[1],".un.hetcount",sep="")) 
hz_hm_s100 <- read.table(paste("s100/",args[1],".hm.hetcount",sep=""))
hz_un_s50 <- read.table(paste("s50/",args[1],".un.hetcount",sep=""))
hz_hm_s50 <- read.table(paste("s50/",args[1],".hm.hetcount",sep=""))
hz_un_s20 <- read.table(paste("s20/",args[1],".un.hetcount",sep=""))
hz_hm_s20 <- read.table(paste("s20/",args[1],".hm.hetcount",sep=""))

#######################################################################################################################################
############################# Plotting repeat content PSMC ############################################################################
#######################################################################################################################################

pdf(paste(args[1],"_repeat_effect_s100.pdf",sep=""),height=10,width=20)
plot(pt_un_s100$V1,(pt_un_s100$V2*10000),log="x",ylim=c(0,60000),type="n", xlab="Years ago", ylab="effective population size",main="Effec
t of repeat on PSMC (s100)")
lines(pt_hm_s100$V1,(pt_hm_s100$V2*10000),type="s",col="blue",lwd=2)
lines(pt_un_s100$V1,(pt_un_s100$V2*10000),type="s",col="orangered",lwd=2)
lines(msmc_hm$left_time_boundary/mu*gen,(1/msmc_hm$lambda_00)/(2*mu),type="s",col="skyblue",lwd=2)
lines(msmc_un$left_time_boundary/mu*gen,(1/msmc_un$lambda_00)/(2*mu),type="s",col="red",lwd=2)
points(pt_un_s100$V1,(pt_un_s100$V2*10000),col="black",pch=20,cex=c(hz_un_s100$V4/10))
legend("topleft",legend=c("PSMC Masked","PSMC Unmasked","MSMC Unmasked","MSMC Masked"),col=c("blue","orangered","red","skyblue"),lwd=3)
dev.off()

pdf(paste(args[1],"_repeat_effect_s50.pdf",sep=""),height=10,width=20)
plot(pt_un_s50$V1,(pt_un_s50$V2*10000),log="x",ylim=c(0,60000),type="n", xlab="Years ago", ylab="effective population size",main="Effect 
of repeat on PSMC (s50)")
lines(pt_hm_s50$V1,(pt_hm_s50$V2*10000),type="s",col="blue",lwd=2)
lines(pt_un_s50$V1,(pt_un_s50$V2*10000),type="s",col="orangered",lwd=2)
lines(msmc_hm$left_time_boundary/mu*gen,(1/msmc_hm$lambda_00)/(2*mu),type="s",col="skyblue",lwd=2)
lines(msmc_un$left_time_boundary/mu*gen,(1/msmc_un$lambda_00)/(2*mu),type="s",col="red",lwd=2)
points(pt_un_s50$V1,(pt_un_s50$V2*10000),col="black",pch=20,cex=c(hz_un_s50$V4/10))
legend("topleft",legend=c("PSMC Masked","PSMC Unmasked","MSMC Unmasked","MSMC Masked"),col=c("blue","orangered","red","skyblue"),lwd=3)
dev.off()

pdf(paste(args[1],"_repeat_effect_s20.pdf",sep=""),height=10,width=20)
plot(pt_un_s20$V1,(pt_un_s20$V2*10000),log="x",ylim=c(0,60000),type="n", xlab="Years ago", ylab="effective population size",main="Effect of repeat on PSMC (s20)")
lines(pt_hm_s20$V1,(pt_hm_s20$V2*10000),type="s",col="blue",lwd=2)
lines(pt_un_s20$V1,(pt_un_s20$V2*10000),type="s",col="orangered",lwd=2)
lines(msmc_hm$left_time_boundary/mu*gen,(1/msmc_hm$lambda_00)/(2*mu),type="s",col="skyblue",lwd=2)
lines(msmc_un$left_time_boundary/mu*gen,(1/msmc_un$lambda_00)/(2*mu),type="s",col="red",lwd=2)
points(pt_un_s20$V1,(pt_un_s20$V2*10000),col="black",pch=20,cex=c(hz_un_s20$V4/10))
legend("topleft",legend=c("PSMC Masked","PSMC Unmasked","MSMC Unmasked","MSMC Masked"),col=c("blue","orangered","red","skyblue"),lwd=2)
dev.off()

pdf(paste(args[1],"_repeat_effect_all.pdf",sep=""),height=10,width=20)
plot(pt_un_s100$V1,(pt_un_s100$V2*10000),log="x",ylim=c(0,60000),type="n", xlab="Years ago", ylab="effective population size",main="Effect of repeat on PSMC (s100)")
lines(pt_hm_s100$V1,(pt_hm_s100$V2*10000),type="s",col="blue",lwd=2)
lines(pt_un_s100$V1,(pt_un_s100$V2*10000),type="s",col="orangered",lwd=2)
lines(msmc_hm$left_time_boundary/mu*gen,(1/msmc_hm$lambda_00)/(2*mu),type="s",col="skyblue",lwd=2)
lines(msmc_un$left_time_boundary/mu*gen,(1/msmc_un$lambda_00)/(2*mu),type="s",col="red",lwd=2)
points(pt_un_s100$V1,(pt_un_s100$V2*10000),col="black",pch=20,cex=c(hz_un_s100$V4/10))
legend("topleft",legend=c("PSMC Masked","PSMC Unmasked","MSMC Unmasked","MSMC Masked"),col=c("blue","orangered","red","skyblue"),lwd=2)
plot(pt_un_s50$V1,(pt_un_s50$V2*10000),log="x",ylim=c(0,60000),type="n", xlab="Years ago", ylab="effective population size",main="Effect of repeat on PSMC (s50)")
lines(pt_hm_s50$V1,(pt_hm_s50$V2*10000),type="s",col="blue",lwd=2)
lines(pt_un_s50$V1,(pt_un_s50$V2*10000),type="s",col="orangered",lwd=2)
lines(msmc_hm$left_time_boundary/mu*gen,(1/msmc_hm$lambda_00)/(2*mu),type="s",col="skyblue",lwd=2)
lines(msmc_un$left_time_boundary/mu*gen,(1/msmc_un$lambda_00)/(2*mu),type="s",col="red",lwd=2)
points(pt_un_s50$V1,(pt_un_s50$V2*10000),col="black",pch=20,cex=c(hz_un_s100$V4/10))
legend("topleft",legend=c("PSMC Masked","PSMC Unmasked","MSMC Unmasked","MSMC Masked"),col=c("blue","orangered","red","skyblue"),lwd=2)
plot(pt_un_s20$V1,(pt_un_s20$V2*10000),log="x",ylim=c(0,60000),type="n", xlab="Years ago", ylab="effective population size",main="Effect of repeat on PSMC (s20)")
lines(pt_hm_s20$V1,(pt_hm_s20$V2*10000),type="s",col="blue",lwd=2)
lines(pt_un_s20$V1,(pt_un_s20$V2*10000),type="s",col="orangered",lwd=2)
lines(msmc_hm$left_time_boundary/mu*gen,(1/msmc_hm$lambda_00)/(2*mu),type="s",col="skyblue",lwd=2)
lines(msmc_un$left_time_boundary/mu*gen,(1/msmc_un$lambda_00)/(2*mu),type="s",col="red",lwd=2)
points(pt_un_s20$V1,(pt_un_s20$V2*10000),col="black",pch=20,cex=c(hz_un_s100$V4/10))
legend("topleft",legend=c("PSMC Masked","PSMC Unmasked","MSMC Unmasked","MSMC Masked"),col=c("blue","orangered","red","skyblue"),lwd=2)
dev.off()

#######################################################################################################################################
##################### Plotting Regression plot between difference in Ne and Repeat content ############################################
#######################################################################################################################################

s100_ne_diff <- pt_un_s100$V2-pt_hm_s100$V2
s50_ne_diff <- pt_un_s50$V2-pt_hm_s50$V2
s20_ne_diff <- pt_un_s20$V2-pt_hm_s20$V2
reg_s100 <- lm(hz_un_s100$V4~s100_ne_diff)
reg_s50 <- lm(hz_un_s50$V4~s50_ne_diff)
reg_s20 <- lm(hz_un_s20$V4~s20_ne_diff)
coeff_s100=coefficients(reg_s100)
coeff_s50=coefficients(reg_s50)
coeff_s20=coefficients(reg_s20)
eq_s100 = paste0("y = ", round(coeff_s100[2],1), " x ", round(coeff_s100[1],1))
eq_s50 = paste0("y = ", round(coeff_s50[2],1), " x ", round(coeff_s50[1],1))
eq_s20 = paste0("y = ", round(coeff_s20[2],1), " x ", round(coeff_s20[1],1))
ken_s100 <- cor.test(hz_un_s100$V4,s100_ne_diff, method="kendall")
ken_s50 <- cor.test(hz_un_s50$V4,s50_ne_diff, method="kendall")
ken_s20 <- cor.test(hz_un_s20$V4,s20_ne_diff, method="kendall")
tau_s100 = paste0("Tau = ", ken_s100$estimate)
tau_s50 = paste0("Tau = ", ken_s50Q$estimate)
tau_s20 = paste0("Tau = ", ken_s20Q$estimate)
pval_s100 = paste0("P-value = ", ken_s100$p.value)
pval_s50 = paste0("P-value = ", ken_s50$p.value)
pval_s20 = paste0("P-value = ", ken_s20$p.value)

pdf(paste(args[1],"_correlation_repeat_s100.pdf",sep=""))
plot(hz_un_s100$V4~s100_ne_diff,main="s100_repeat_content_correlation",ylab="Repeat content",xlab="Difference in Ne",pch=20)
abline(reg_s100, col="blue",lwd=2,lty=4)
text(mean(s100_ne_diff)+0.1, 0.2 * max(hz_un_s100$V4) + 0, labels=eq_s100, srt=0.2, col = "red")
text(mean(s100_ne_diff)+0.1,max(hz_un_s100$V4) + 0, labels=tau_s100, col = "red")
text(mean(s100_ne_diff)+0.1,max(hz_un_s100$V4) - 4, labels=pval_s100, col = "red")
dev.off()

pdf(paste(args[1],"_correlation_repeat_s50.pdf",sep=""))
plot(hz_un_s50$V4~s50_ne_diff,main="s50_repeat_content_correlation",ylab="Repeat content",xlab="Difference in Ne",pch=20)
abline(reg_s50, col="blue",lwd=2,lty=4)
text(mean(s50_ne_diff)+0.1, 0.2 * max(hz_un_s50$V4) + 0, labels=eq_s50, srt=0.2, col = "red")
text(mean(s50_ne_diff)+0.1,max(hz_un_s50$V4) + 0, labels=tau_s50, col = "red")
text(mean(s50_ne_diff)+0.1,max(hz_un_s50$V4) - 4, labels=pval_s50, col = "red")
dev.off()

pdf(paste(args[1],"_correlation_repeat_s20.pdf",sep=""))
plot(hz_un_s20$V4~s20_ne_diff,main="s20_repeat_content_correlation",ylab="Repeat content",xlab="Difference in Ne",pch=20)
abline(reg_s20, col="blue",lwd=2,lty=4)
text(mean(s20_ne_diff)+0.1, 0.2 * max(hz_un_s20$V4) + 0, labels=eq_s20, srt=0.2, col = "red")
text(mean(s20_ne_diff)+0.1,max(hz_un_s20$V4) + 0, labels=tau_s20, col = "red")
text(mean(s20_ne_diff)+0.1,max(hz_un_s20$V4) - 4, labels=pval_s20, col = "red")
dev.off()

pdf(paste(args[1],"_correlation_repeat_all.pdf",sep=""))
plot(hz_un_s100$V4~s100_ne_diff,main="s100_repeat_content_correlation",ylab="Repeat content",xlab="Difference in Ne",pch=20)
abline(reg_s100, col="blue",lwd=2,lty=4)
text(mean(s100_ne_diff)+0.1, 0.2 * max(hz_un_s100$V4) + 0, labels=eq_s100, srt=0.2, col = "red")
text(mean(s100_ne_diff)+0.1,max(hz_un_s100$V4) + 0, labels=tau_s100, col = "red")
text(mean(s100_ne_diff)+0.1,max(hz_un_s100$V4) - 4, labels=pval_s100, col = "red")
plot(hz_un_s50$V4~s50_ne_diff,main="s50_repeat_content_correlation",ylab="Repeat content",xlab="Difference in Ne",pch=20)
abline(reg_s50, col="blue",lwd=2,lty=4)
text(mean(s50_ne_diff)+0.1, 0.2 * max(hz_un_s50$V4) + 0, labels=eq_s50, srt=0.2, col = "red")
text(mean(s50_ne_diff)+0.1,max(hz_un_s50$V4) + 0, labels=tau_s50, col = "red")
text(mean(s50_ne_diff)+0.1,max(hz_un_s50$V4) - 4, labels=pval_s50, col = "red")
plot(hz_un_s20$V4~s20_ne_diff,main="s20_repeat_content_correlation",ylab="Repeat content",xlab="Difference in Ne",pch=20)
abline(reg_s20, col="blue",lwd=2,lty=4)
text(mean(s20_ne_diff)+0.1, 0.2 * max(hz_un_s20$V4) + 0, labels=eq_s20, srt=0.2, col = "red")
text(mean(s20_ne_diff)+0.1,max(hz_un_s20$V4) + 0, labels=tau_s20, col = "red")
text(mean(s20_ne_diff)+0.1,max(hz_un_s20$V4) - 4, labels=pval_s20, col = "red")
dev.off()

#######################################################################################################################################

data.frame2matrix = function(data, rowtitle, coltitle, datatitle, 
                             rowdecreasing = FALSE, coldecreasing = FALSE,
                             default_value = NA) {

  # check, whether titles exist as columns names in the data.frame data
  if ( (!(rowtitle%in%names(data))) 
       || (!(coltitle%in%names(data))) 
       || (!(datatitle%in%names(data))) ) {
    stop('data.frame2matrix: bad row-, col-, or datatitle.')
  }

  # get number of rows in data
  ndata = dim(data)[1]

  # extract rownames and colnames for the matrix from the data.frame
  rownames = sort(unique(data[[rowtitle]]), decreasing = rowdecreasing)
  nrows = length(rownames)
  colnames = sort(unique(data[[coltitle]]), decreasing = coldecreasing)
  ncols = length(colnames)

  # initialize the matrix
  out_matrix = matrix(NA, 
                      nrow = nrows, ncol = ncols,
                      dimnames=list(rownames, colnames))

  # iterate rows of data
  for (i1 in 1:ndata) {
    # get matrix-row and matrix-column indices for the current data-row
    iR = which(rownames==data[[rowtitle]][i1])
    iC = which(colnames==data[[coltitle]][i1])

    # throw an error if the matrix entry (iR,iC) is already filled.
    if (!is.na(out_matrix[iR, iC])) stop('data.frame2matrix: double entry in data.frame')
    out_matrix[iR, iC] = data[[datatitle]][i1]
  }

  # set empty matrix entries to the default value
  out_matrix[is.na(out_matrix)] = default_value

  # return matrix
  return(out_matrix)

}

########################################################################################################################################


#########################################################################################################################################
####################################### INTERVAL WISE REPEAT DISTRIBUTION ###############################################################
#########################################################################################################################################

temp = list.files("repclass",pattern=".rc.bed",full.names=TRUE)
for (i in temp){
read.table(file=i,header=FALSE)-> rc
rc$V3-rc$V2 -> rc$diff
as.data.frame(aggregate(rc$diff,list(rc$V4),sum))-> aggr_rc
aggr_rc$percent_cont <- (aggr_rc$x/sum(aggr_rc$x))*100
data.frame(Percent=aggr_rc$percent_cont,Rep_class=aggr_rc$Group.1,atomInt=rep(strsplit(i, "[.]")[[1]][1],length(aggr_rc$Group.1)))->X
write.table(file="AI_perc.txt",X,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
}
rep_col <- c("orangered","blue","plum1","hotpink","darkolivegreen","peru","olivedrab1","orange2","darkseagreen1","palevioletred1","navy","khaki","gray50","cornflowerblue","yellowgreen","chartreuse4","cadetblue4","burlywood4","brown4","blueviolet","bisque4")
read.table(file="AI_perc.txt",header=FALSE)->M
myMatrix = data.frame2matrix(M, 'V2', 'V3', 'V1')
myMatrix[is.na(myMatrix)] <- 0
pdf("AI_perc.pdf",height=15,width=20)
par(mar = c(8.1,4.1, 3.1, 3.1),xpd=TRUE)
barplot(myMatrix, col=rep_col, main="Percent of intervals covered by repeat classes",xlab="Atomic Interval",ylab="Abundance")
legend("bottom", inset=c(0,-0.31),legend=c(row.names(myMatrix)),fill=rep_col,ncol = 5, border=NA,box.lty=2,cex=0.6)
dev.off()

#########################################################################################################################################
################################# GENOMEWIDE REPEAT DISTRIBUTION ACROSS INTERVALS #######################################################
#########################################################################################################################################

temp = list.files("repclass",pattern=".rc.bed",full.names=TRUE)
for (i in temp){
read.table(file=i,header=FALSE)-> rc
rc$V3-rc$V2 -> rc$diff
as.data.frame(aggregate(rc$diff,list(rc$V4),sum))-> aggr_rc
aggr_rc$percent_geno <- paste(aggr_rc$x/args[4])*100
data.frame(Percent=aggr_rc$percent_geno,Rep_class=aggr_rc$Group.1,atomInt=rep(strsplit(i, "[.]")[[1]][1],length(aggr_rc$Group.1)))->Y
write.table(file="AI_perc_geno.txt",Y,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
}
read.table(file="AI_perc_geno.txt",header=FALSE)-> R
genMatrix = data.frame2matrix(R, 'V2', 'V3', 'V1')
genMatrix[is.na(genMatrix)] <- 0
pdf("AI_perc_geno.pdf",height=15,width=20)
par(mar = c(8.1,4.1, 3.1, 3.1),xpd=TRUE)
barplot(genMatrix, col=MyCol, main="Percent of genome covered by intervals and repeat classes",xlab="Atomic Interval",ylab="genome percent")
legend("bottom", inset=c(0,-0.31),legend=c(row.names(genMatrix)),fill=rep_col , ncol = 5, border=NA, box.lty=2, cex=0.6)
dev.off()

########################################################################################################################################
#################################### Plotting Repeat Class Ts/Tv and Hety Plots #######################################################
#######################################################################################################################################

O_s20 <- read.table(paste("s20/",args[1],".rep.profile",sep=""),header=TRUE)
O_s50 <- read.table(paste("s50/",args[1],".rep.profile",sep=""),header=TRUE)
O_s100 <- read.table(paste("s100/",args[1],".rep.profile",sep=""),header=TRUE)
N_s20 <- O_s20[is.finite(O_s20$Ts.Tv_ratio),]
N_s20 <- N_s20[is.finite(N_s20$Hety_percent),]
M_s20 <- N_s20[order(N_s20$Atomic_Interval),]
N_s50 <- O_s50[is.finite(O_s50$Ts.Tv_ratio),]
N_s50 <- N_s50[is.finite(N_s50$Hety_percent),]
M_s50 <- N_s50[order(N_s50$Atomic_Interval),]
N_s100 <- O_s100[is.finite(O_s100$Ts.Tv_ratio),]
N_s100 <- N_s100[is.finite(N_s100$Hety_percent),]
M_s100 <- N_s100[order(N_s100$Atomic_Interval),]

pdf(args[1],".rep_tstv.s100.pdf",height=15,width=20)
par(mfrow=c(4,5))
N <- O_s100[is.finite(O_s100$Ts.Tv_ratio),]
N <- N[is.finite(N$Hety_percent),]
M <- N[order(N$Atomic_Interval),]
colcount<-1
for (repclass in unique(M$Repeatclass)) {
colcount<-colcount+1
par(mar=c(5, 4, 4, 5))
plot(M$Atomic_Interval[M$Repeatclass==repclass],M$Hety_percent[M$Repeatclass==repclass],axes=FALSE,ylim=c(0,2),xlim=c(0,63),type="l",main=repclass,xlab="Atomic Interval",ylab="",col="black")
axis(1, xlim=c(0,63),col="black",las=1)
axis(2, ylim=c(0,2),col="black",las=1)
mtext("hety %",side=2,line=2.5)
box()
par(new=TRUE)
plot(M$Atomic_Interval[M$Repeatclass==repclass],M$Ts.Tv_ratio[M$Repeatclass==repclass],axes=FALSE,ylim=c(0,32),type="l",xlab="",ylab="",col="red")
mtext("Ts/Tv",side=4,col="red",line=2.5)
axis(4, ylim=c(0,32), col="red",col.axis="red",las=1)
}
dev.off()

pdf(args[1],".rep_tstv.s50.pdf",height=15,width=20)
par(mfrow=c(4,5))
N <- O_s50[is.finite(O_s50$Ts.Tv_ratio),]
N <- N[is.finite(N$Hety_percent),]
M <- N[order(N$Atomic_Interval),]
colcount<-1
for (repclass in unique(M$Repeatclass)) {
colcount<-colcount+1
par(mar=c(5, 4, 4, 5))
plot(M$Atomic_Interval[M$Repeatclass==repclass],M$Hety_percent[M$Repeatclass==repclass],axes=FALSE,ylim=c(0,2),xlim=c(0,63),type="l",main=repclass,xlab="Atomic Interval",ylab="",col="black")
axis(1, xlim=c(0,63),col="black",las=1)
axis(2, ylim=c(0,2),col="black",las=1)
mtext("hety %",side=2,line=2.5)
box()
par(new=TRUE)
plot(M$Atomic_Interval[M$Repeatclass==repclass],M$Ts.Tv_ratio[M$Repeatclass==repclass],axes=FALSE,ylim=c(0,32),type="l",xlab="",ylab="",col="red")
mtext("Ts/Tv",side=4,col="red",line=2.5)
axis(4, ylim=c(0,32), col="red",col.axis="red",las=1)
}
dev.off()

pdf(args[1],".rep_tstv.s20.pdf",height=15,width=20)
par(mfrow=c(4,5))
N <- O_s20[is.finite(O_s20$Ts.Tv_ratio),]
N <- N[is.finite(N$Hety_percent),]
M <- N[order(N$Atomic_Interval),]
colcount<-1
for (repclass in unique(M$Repeatclass)) {
colcount<-colcount+1
par(mar=c(5, 4, 4, 5))
plot(M$Atomic_Interval[M$Repeatclass==repclass],M$Hety_percent[M$Repeatclass==repclass],axes=FALSE,ylim=c(0,2),xlim=c(0,63),type="l",main=repclass,xlab="Atomic Interval",ylab="",col="black")
axis(1, xlim=c(0,63),col="black",las=1)
axis(2, ylim=c(0,2),col="black",las=1)
mtext("hety %",side=2,line=2.5)
box()
par(new=TRUE)
plot(M$Atomic_Interval[M$Repeatclass==repclass],M$Ts.Tv_ratio[M$Repeatclass==repclass],axes=FALSE,ylim=c(0,32),type="l",xlab="",ylab="",col="red")
mtext("Ts/Tv",side=4,col="red",line=2.5)
axis(4, ylim=c(0,32), col="red",col.axis="red",las=1)
}
dev.off()

#######################################################################################################################################
############################################ Plotting repeat class PSMC ###############################################################
#######################################################################################################################################

rep_col <- c("plum1","hotpink","darkolivegreen","peru","olivedrab1","orange2","darkseagreen1","palevioletred1","navy","khaki","gray50","cornflowerblue","yellowgreen","chartreuse4","cadetblue4","burlywood4","brown4","blueviolet","bisque4")
colcount <- 0
rep = list.files("repclass",pattern=".0.txt",full.names=TRUE)
repclass <- sapply(strsplit(sapply(strsplit(rep,"[.]"), `[`, 1),"[/]"), `[`, 2)
pdf(paste(args[1],"_repclass_PSMC_s100.pdf",sep=""),height=10,width=20)
plot(pt_un_s100$V1,(pt_un_s100$V2*10000),log="x",ylim=c(0,60000),type="n", xlab="Years ago", ylab="effective population size",main="Effect of repeat on PSMC (s100)")
lines(pt_un_s100$V1,(pt_un_s100$V2*10000),type="s",col="orangered",lwd=2)
lines(pt_hm_s100$V1,(pt_hm_s100$V2*10000),type="s",col="blue",lwd=2)
for (i in rep){
colcount <- colcount+1
read.table(file=i,header=FALSE)-> repc
lines(repc$V1,(repc$V2*10000),type="s",col=rep_col[colcount],lwd=2)
}
legend("top",legend=repclass,fill=rep_col,ncol = 5,cex=1.75)

#######################################################################################################################################
################################### Plotting No. of recombinations across Atomic Intervals ############################################
#######################################################################################################################################

nrec <- list.files(pattern=".nrec")
for (i in nrec){
read.table(file=i,header=FALSE)-> M
length(M$V1)/64->itercount
ai<-64
t(matrix(M$V1,nrow = ai,ncol = itercount))->N
pdf(paste(i,"_RecQC.pdf",sep=""),width=15,height=20)
par(mfrow=c(4,4))
for(itK in c(1:ai)){
barplot(N[,itK],beside=TRUE,main=itK,ylab="No.of recombinations",xlab="No .of iteration",font.lab=2,cex.lab=1.2,cex.main=1.5)
abline(h=10,v=20,col="red", lty =4,lwd=2)
box()
}
dev.off()
}

#######################################################################################################################################
############################### Plotting Lengths of sequences across Atomic Intervals #################################################
#######################################################################################################################################

rlen <- list.files(pattern=".lengths")
pdf("RepLength.pdf",width=15,height=20)
par(mfrow=c(2,2))
for (i in rlen){
read.table(file=i,header=FALSE)-> M
boxplot(M$V2~M$V1,log="y",main=i,xlab="Atomic Intervals",ylab="log(Lengths)",font.lab=2,cex.lab=1.2,cex.main=1.5)
box()
}
dev.off()

#######################################################################################################################################



