
setwd("/Users/lauradonohue/Hi-ChIP/4C/05162017_A431-1-CTRLi_p63i-igv")
library(Hmisc)
library(zoo)

######################################
######################################

# Read in the data

center_coord=106530000
resolution=10000
table1=read.table("4C_PRDM1_igv_A431-1-CTRLi_10000.matrix.bed.bed")
table2=read.table("4C_PRDM1_igv_A431-1-p63i_10000.matrix.bed.bed")
#table3=read.table("H9_merge_chr17_70110000_10kb.matrix.bed")
#table4=read.table("WT33_r1_chr17_70110000_10kb.matrix.bed")


######################################
######################################

# Process data

#Normalization factors: take the number of read mapped for whole library (HiC-Pro->data->sample-> *.mergestat->valid interaction rmdup ->NUMBER and divide by million)

# Normalization factors
table1_norm=155.074633 # A431-1-CTRLi HiChIP
table2_norm=117.114388  # A431-1-p63i HiChIP
#table3_norm=155.074633 #5 CTRLi HiChIP
#table4_norm=33.187307 #1m


# Alternate norm factors

table1_norm = sum(table1[,7])/1e3
table2_norm = sum(table2[,7])/1e3


################


# Pick your smoothing window (1 means no smoothing, 2 averages over 2 adjacent bins, etc. )
rollmean_win =3

# Process dataset 1
counts1_left=table1[table1[,2]== center_coord,c(5,7)]
counts1_right=table1[table1[,5]== center_coord,c(2,7)]
colnames(counts1_left) = c(1,2)
colnames(counts1_right) = c(1,2)
counts1_merge=rbind(counts1_left, counts1_right)
counts1_merge_ord = counts1_merge[order(counts1_merge[,1]),]
counts1_merge_ord[,2] = counts1_merge_ord[,2]/table1_norm
counts1_merge_ord[,1] = counts1_merge_ord[,1] + resolution/2
line1 = rollmean(counts1_merge_ord[,2], rollmean_win,fill=c(0,0,0))


# Process dataset 2
counts2_left=table2[table2[,2]== center_coord,c(5,7)]
counts2_right=table2[table2[,5]== center_coord,c(2,7)]
colnames(counts2_left) = c(1,2)
colnames(counts2_right) = c(1,2)
counts2_merge=rbind(counts2_left, counts2_right)
counts2_merge_ord = counts2_merge[order(counts2_merge[,1]),]
counts2_merge_ord[,2] = counts2_merge_ord[,2]/table2_norm
counts2_merge_ord[,1] = counts2_merge_ord[,1] + resolution/2
line2 = rollmean(counts2_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset3
counts3_left=table3[table3[,2]== center_coord,c(5,7)]
counts3_right=table3[table3[,5]== center_coord,c(2,7)]
colnames(counts3_left) = c(1,2)
colnames(counts3_right) = c(1,2)
counts3_merge=rbind(counts3_left, counts3_right)
counts3_merge_ord = counts3_merge[order(counts3_merge[,1]),]
counts3_merge_ord[,2] = counts3_merge_ord[,2]/table3_norm
counts3_merge_ord[,1] = counts3_merge_ord[,1] + resolution/2
line3 = rollmean(counts3_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset4
counts4_left=table4[table4[,2]== center_coord,c(5,7)]
counts4_right=table4[table4[,5]== center_coord,c(2,7)]
colnames(counts4_left) = c(1,2)
colnames(counts4_right) = c(1,2)
counts4_merge=rbind(counts4_left, counts4_right)
counts4_merge_ord = counts4_merge[order(counts4_merge[,1]),]
counts4_merge_ord[,2] = counts4_merge_ord[,2]/table4_norm
counts4_merge_ord[,1] = counts4_merge_ord[,1] + resolution/2
line4 = rollmean(counts4_merge_ord[,2], rollmean_win,fill=c(0,0,0))



######################################
######################################

# Plot

xlim_minus=1e6 # bp upstream of anchor
xlim_plus=1e6 # bp downstream of anchor

ylim_par=c(-0.01,2) # y axis limits

lwd_par = 2 # line thickness
xlab="Position on chr6"
ylab="Reads/1e6 filtered reads"
plot_name=paste("Virtual 4C anchor 10kb bin ",center_coord,sep="")
xlim_par= c((center_coord)-xlim_minus,(center_coord)+xlim_plus)

plot(counts1_merge_ord[,1],line1,type="l",xlab= xlab,ylab=ylab,lwd= lwd_par,xlim=xlim_par,main=plot_name,ylim= ylim_par,cex.main=0.7,col="blue",xaxs="i")
#points(counts1_merge_ord[,1], counts1_merge_ord[,2],pch=19,cex=0.4,col="blue")

lines(counts2_merge_ord[,1],line2,type="l",col="red",lwd=lwd_par)
#points(counts2_merge_ord[,1], counts2_merge_ord[,2],pch=19,col="blue",lwd=lwd_par,cex=0.4)

#lines(counts3_merge_ord[,1],line3,type="l",col="purple",lwd=lwd_par)
#points(counts3_merge_ord[,1], counts3_merge_ord[,2],pch=19,col="red",lwd=lwd_par,cex=0.4)

#lines(counts4_merge_ord[,1],line4,type="l",col="darkgreen",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)


abline(v= center_coord + resolution/2,lty=2,col="gray")
#legend("topright",c("0818","1209","H9","WT33"), col=c("blue","red","purple","darkgreen"),inset=0.03,lwd=lwd_par,cex=0.7)

legend("topright",c("A431-1-CTRLi", "A431-1-p63i"), col=c("blue", "red"),inset=0.03,lwd=lwd_par,cex=0.7)
