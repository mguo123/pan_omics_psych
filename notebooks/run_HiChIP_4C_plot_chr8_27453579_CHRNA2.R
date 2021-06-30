# setwd("/Users/lauradonohue/Khavari_Bioinformatics/Project/Wernig_collab/HiChIP/4C_plots")
# setwd("/Users/lauradonohue/Khavari_Bioinformatics/Project/Wernig_collab/HiChIP/4C_plots")
library(Hmisc)
library(zoo)

######################################
######################################
# read in dataframes and set save path

# set variables
samples = c("Astro_B2","H9D10-B1","SL-B2","SLC-B1")


# Read in the data
gene_sym="CHRNA2"
chrom=8
center_coord=27453579
resolution=5000
save_file = paste0(save_dir,gene_sym, '_4tissues.pdf')
save_file
# samples = c("Astro_B2","H9-B2","H9D28-B1","SLC-B1")
table1=read.table(paste0(dir_4C,"4C_", gene_sym, '_', samples[1],"_5000.matrix.bed"))
table2=read.table(paste0(dir_4C,"4C_", gene_sym, '_', samples[2],"_5000.matrix.bed"))



center_coord=27450000
resolution=5000
samples = c('Astro_B1', 'Astro_B2', 'H9-B1','H9-B2','H9D10-B1','H9D10-B1',
            'H9D10-B1','H9D10-B1',)
table1=read.table("4C_chr8_27453579_CHRNA2_Astro_B1_5000.matrix.bed")
table2=read.table("4C_chr8_27453579_CHRNA2_Astro_B2_5000.matrix.bed")
table3=read.table("4C_chr8_27453579_CHRNA2_H9-B1_5000.matrix.bed")
table4=read.table("4C_chr8_27453579_CHRNA2_H9-B2_5000.matrix.bed")
table5=read.table("4C_chr8_27453579_CHRNA2_H9D10-B1_5000.matrix.bed")
table6=read.table("4C_chr8_27453579_CHRNA2_H9D10-B2_5000.matrix.bed")
table7=read.table("4C_chr8_27453579_CHRNA2_H9D28-B1_5000.matrix.bed")
table8=read.table("4C_chr8_27453579_CHRNA2_H9D28-B2_5000.matrix.bed")
table9=read.table("4C_chr8_27453579_CHRNA2_SL-B1_5000.matrix.bed")
table10=read.table("4C_chr8_27453579_CHRNA2_SL-B2_5000.matrix.bed")
table11=read.table("4C_chr8_27453579_CHRNA2_SLC-B1_5000.matrix.bed")
table12=read.table("4C_chr8_27453579_CHRNA2_SLC-B2_5000.matrix.bed")
table13=read.table("4C_chr8_27453579_CHRNA2_H9D2-B1_5000.matrix.bed")
table14=read.table("4C_chr8_27453579_CHRNA2_H9D2-B2_5000.matrix.bed")
table15=read.table("4C_chr8_27453579_CHRNA2_SLD2-B1_5000.matrix.bed")
table16=read.table("4C_chr8_27453579_CHRNA2_SLD2-B2_5000.matrix.bed")
table17=read.table("4C_chr8_27453579_CHRNA2_SLCD2-B1_5000.matrix.bed")
table18=read.table("4C_chr8_27453579_CHRNA2_SLCD2-B2_5000.matrix.bed")


######################################
######################################

# Process data

# Normalization factors: take the number of read mapped for whole library (HiC-Pro->data->sample-> *.mergestat->valid interaction rmdup ->NUMBER and divide by million)

table1_norm=224.66504
table2_norm=209.890877
table3_norm=141.426659
table4_norm=116.35082
table5_norm=95.190714
table6_norm=88.661039
table7_norm=148.778866
table8_norm=135.892507
table9_norm=180.870221
table10_norm=272.485208
table11_norm=135.989382
table12_norm=162.691811
table13_norm=340.200815
table14_norm=348.373024
table15_norm=319.537418
table16_norm=318.400529
table17_norm=308.753333
table18_norm=201.431031

# table1_norm=1 
# table2_norm = 1
# table3_norm = 1
# table4_norm = 1
# table5_norm = 1
# table6_norm = 1
# table7_norm = 1
# table8_norm = 1
# table9_norm = 1
# table10_norm = 1
# table11_norm = 1
# table12_norm = 1
# table13_norm = 1
# table14_norm = 1
# table15_norm = 1
# table16_norm = 1
# table17_norm = 1


######################################
######################################

# Alternate norm factors

#table1_norm = sum(table1[,7])/1e3
#table2_norm = sum(table2[,7])/1e3
#table3_norm = sum(table3[,7])/1e3
#table4_norm = sum(table4[,7])/1e3

################


# Pick your smoothing window (1 means no smoothing, 2 averages over 2 adjacent bins, etc. )
rollmean_win =2

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

# Process dataset 5
counts5_left=table5[table5[,2]== center_coord,c(5,7)]
counts5_right=table5[table5[,5]== center_coord,c(2,7)]
colnames(counts5_left) = c(1,2)
colnames(counts5_right) = c(1,2)
counts5_merge=rbind(counts5_left, counts5_right)
counts5_merge_ord = counts5_merge[order(counts5_merge[,1]),]
counts5_merge_ord[,2] = counts5_merge_ord[,2]/table5_norm
counts5_merge_ord[,1] = counts5_merge_ord[,1] + resolution/2
line5 = rollmean(counts5_merge_ord[,2], rollmean_win,fill=c(0,0,0))


# Process dataset 6
counts6_left=table6[table6[,2]== center_coord,c(5,7)]
counts6_right=table6[table6[,5]== center_coord,c(2,7)]
colnames(counts6_left) = c(1,2)
colnames(counts6_right) = c(1,2)
counts6_merge=rbind(counts6_left, counts6_right)
counts6_merge_ord = counts6_merge[order(counts6_merge[,1]),]
counts6_merge_ord[,2] = counts6_merge_ord[,2]/table6_norm
counts6_merge_ord[,1] = counts6_merge_ord[,1] + resolution/2
line6 = rollmean(counts6_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset7
counts7_left=table7[table7[,2]== center_coord,c(5,7)]
counts7_right=table7[table7[,5]== center_coord,c(2,7)]
colnames(counts7_left) = c(1,2)
colnames(counts7_right) = c(1,2)
counts7_merge=rbind(counts7_left, counts7_right)
counts7_merge_ord = counts7_merge[order(counts7_merge[,1]),]
counts7_merge_ord[,2] = counts7_merge_ord[,2]/table7_norm
counts7_merge_ord[,1] = counts7_merge_ord[,1] + resolution/2
line7 = rollmean(counts7_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset8
counts8_left=table8[table8[,2]== center_coord,c(5,7)]
counts8_right=table8[table8[,5]== center_coord,c(2,7)]
colnames(counts8_left) = c(1,2)
colnames(counts8_right) = c(1,2)
counts8_merge=rbind(counts8_left, counts8_right)
counts8_merge_ord = counts8_merge[order(counts8_merge[,1]),]
counts8_merge_ord[,2] = counts8_merge_ord[,2]/table8_norm
counts8_merge_ord[,1] = counts8_merge_ord[,1] + resolution/2
line8 = rollmean(counts8_merge_ord[,2], rollmean_win,fill=c(0,0,0))

# Process dataset 9
counts9_left=table9[table9[,2]== center_coord,c(5,7)]
counts9_right=table9[table1[,5]== center_coord,c(2,7)]
colnames(counts9_left) = c(1,2)
colnames(counts9_right) = c(1,2)
counts9_merge=rbind(counts9_left, counts9_right)
counts9_merge_ord = counts9_merge[order(counts9_merge[,1]),]
counts9_merge_ord[,2] = counts9_merge_ord[,2]/table9_norm
counts9_merge_ord[,1] = counts9_merge_ord[,1] + resolution/2
line9 = rollmean(counts9_merge_ord[,2], rollmean_win,fill=c(0,0,0))


# Process dataset 10
counts10_left=table10[table10[,2]== center_coord,c(5,7)]
counts10_right=table10[table10[,5]== center_coord,c(2,7)]
colnames(counts10_left) = c(1,2)
colnames(counts10_right) = c(1,2)
counts10_merge=rbind(counts10_left, counts10_right)
counts10_merge_ord = counts10_merge[order(counts10_merge[,1]),]
counts10_merge_ord[,2] = counts10_merge_ord[,2]/table10_norm
counts10_merge_ord[,1] = counts10_merge_ord[,1] + resolution/2
line10 = rollmean(counts10_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset11
counts11_left=table11[table11[,2]== center_coord,c(5,7)]
counts11_right=table11[table11[,5]== center_coord,c(2,7)]
colnames(counts11_left) = c(1,2)
colnames(counts11_right) = c(1,2)
counts11_merge=rbind(counts11_left, counts11_right)
counts11_merge_ord = counts11_merge[order(counts11_merge[,1]),]
counts11_merge_ord[,2] = counts11_merge_ord[,2]/table11_norm
counts11_merge_ord[,1] = counts11_merge_ord[,1] + resolution/2
line11 = rollmean(counts11_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset12
counts12_left=table12[table12[,2]== center_coord,c(5,7)]
counts12_right=table12[table12[,5]== center_coord,c(2,7)]
colnames(counts12_left) = c(1,2)
colnames(counts12_right) = c(1,2)
counts12_merge=rbind(counts12_left, counts12_right)
counts12_merge_ord = counts12_merge[order(counts12_merge[,1]),]
counts12_merge_ord[,2] = counts12_merge_ord[,2]/table12_norm
counts12_merge_ord[,1] = counts12_merge_ord[,1] + resolution/2
line12 = rollmean(counts12_merge_ord[,2], rollmean_win,fill=c(0,0,0))
# 

# Process dataset 13
counts13_left=table13[table13[,2]== center_coord,c(5,7)]
counts13_right=table13[table13[,5]== center_coord,c(2,7)]
colnames(counts13_left) = c(1,2)
colnames(counts13_right) = c(1,2)
counts13_merge=rbind(counts13_left, counts13_right)
counts13_merge_ord = counts13_merge[order(counts13_merge[,1]),]
counts13_merge_ord[,2] = counts13_merge_ord[,2]/table13_norm
counts13_merge_ord[,1] = counts13_merge_ord[,1] + resolution/2
line13 = rollmean(counts13_merge_ord[,2], rollmean_win,fill=c(0,0,0))


# Process dataset 14
counts14_left=table14[table14[,2]== center_coord,c(5,7)]
counts14_right=table14[table14[,5]== center_coord,c(2,7)]
colnames(counts14_left) = c(1,2)
colnames(counts14_right) = c(1,2)
counts14_merge=rbind(counts14_left, counts14_right)
counts14_merge_ord = counts14_merge[order(counts14_merge[,1]),]
counts14_merge_ord[,2] = counts14_merge_ord[,2]/table14_norm
counts14_merge_ord[,1] = counts14_merge_ord[,1] + resolution/2
line14 = rollmean(counts14_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset15
counts15_left=table15[table15[,2]== center_coord,c(5,7)]
counts15_right=table15[table15[,5]== center_coord,c(2,7)]
colnames(counts15_left) = c(1,2)
colnames(counts15_right) = c(1,2)
counts15_merge=rbind(counts15_left, counts15_right)
counts15_merge_ord = counts15_merge[order(counts15_merge[,1]),]
counts15_merge_ord[,2] = counts15_merge_ord[,2]/table15_norm
counts15_merge_ord[,1] = counts15_merge_ord[,1] + resolution/2
line15 = rollmean(counts15_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset16

counts16_left=table16[table16[,2]== center_coord,c(5,7)]
counts16_right=table16[table16[,5]== center_coord,c(2,7)]
colnames(counts16_left) = c(1,2)
colnames(counts16_right) = c(1,2)
counts16_merge=rbind(counts16_left, counts16_right)
counts16_merge_ord = counts16_merge[order(counts16_merge[,1]),]
counts16_merge_ord[,2] = counts16_merge_ord[,2]/table16_norm
counts16_merge_ord[,1] = counts16_merge_ord[,1] + resolution/2
line16 = rollmean(counts16_merge_ord[,2], rollmean_win,fill=c(0,0,0))

#Dataset17

counts17_left=table17[table17[,2]== center_coord,c(5,7)]
counts17_right=table17[table17[,5]== center_coord,c(2,7)]
colnames(counts17_left) = c(1,2)
colnames(counts17_right) = c(1,2)
counts17_merge=rbind(counts17_left, counts17_right)
counts17_merge_ord = counts17_merge[order(counts17_merge[,1]),]
counts17_merge_ord[,2] = counts17_merge_ord[,2]/table17_norm
counts17_merge_ord[,1] = counts17_merge_ord[,1] + resolution/2
line17 = rollmean(counts17_merge_ord[,2], rollmean_win,fill=c(0,0,0))

# 
#Dataset18

counts18_left=table18[table18[,2]== center_coord,c(5,7)]
counts18_right=table18[table18[,5]== center_coord,c(2,7)]
colnames(counts18_left) = c(1,2)
colnames(counts18_right) = c(1,2)
counts18_merge=rbind(counts18_left, counts18_right)
counts18_merge_ord = counts18_merge[order(counts18_merge[,1]),]
counts18_merge_ord[,2] = counts18_merge_ord[,2]/table18_norm
counts18_merge_ord[,1] = counts18_merge_ord[,1] + resolution/2
line18 = rollmean(counts18_merge_ord[,2], rollmean_win,fill=c(0,0,0))
# 
######################################
######################################

# Plot

#xlim_minus=1.6e5 # bp upstream of anchor
xlim_minus=1.7e5 # bp upstream of anchor
xlim_plus=4.5e4 # bp downstream of anchor

ylim_par=c(-0.01,0.2) # y axis limits

lwd_par = 2 # line thickness
xlab="Position on chr1"
ylab="Filtered reads"
plot_name=paste("Virtual 4C anchor 5kb bin ",center_coord,sep="")
xlim_par= c((center_coord)-xlim_minus,(center_coord)+xlim_plus)

plot(counts1_merge_ord[,1],line1,type="l",xlab= xlab,ylab=ylab,lwd= lwd_par,xlim=xlim_par,main=plot_name,ylim= ylim_par,cex.main=0.7,col="#4E6AE1",xaxs="i")
#points(counts1_merge_ord[,1], counts1_merge_ord[,2],pch=19,cex=0.4,col="blue")

lines(counts2_merge_ord[,1],line2,type="l",col="#4E6AE1",lwd=lwd_par)
#points(counts2_merge_ord[,1], counts2_merge_ord[,2],pch=19,col="blue",lwd=lwd_par,cex=0.4)

lines(counts3_merge_ord[,1],line3,type="l",col="#D5AAE1",lwd=lwd_par)
#points(counts3_merge_ord[,1], counts3_merge_ord[,2],pch=19,col="red",lwd=lwd_par,cex=0.4)

lines(counts4_merge_ord[,1],line4,type="l",col="#D5AAE1",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)

lines(counts5_merge_ord[,1],line5,type="l",col="#740792",lwd=lwd_par)
#points(counts2_merge_ord[,1], counts2_merge_ord[,2],pch=19,col="blue",lwd=lwd_par,cex=0.4)

lines(counts6_merge_ord[,1],line6,type="l",col="#740792",lwd=lwd_par)
#points(counts3_merge_ord[,1], counts3_merge_ord[,2],pch=19,col="red",lwd=lwd_par,cex=0.4)

lines(counts7_merge_ord[,1],line7,type="l",col="#58066F",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)

lines(counts8_merge_ord[,1],line8,type="l",col="#58066F",lwd=lwd_par)
#points(counts2_merge_ord[,1], counts2_merge_ord[,2],pch=19,col="blue",lwd=lwd_par,cex=0.4)

lines(counts9_merge_ord[,1],line9,type="l",col="#FF8CB8",lwd=lwd_par)
#points(counts3_merge_ord[,1], counts3_merge_ord[,2],pch=19,col="red",lwd=lwd_par,cex=0.4)

lines(counts10_merge_ord[,1],line10,type="l",col="#FF8CB8",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)

lines(counts11_merge_ord[,1],line11,type="l",col="#85BF67",lwd=lwd_par)
#points(counts2_merge_ord[,1], counts2_merge_ord[,2],pch=19,col="blue",lwd=lwd_par,cex=0.4)

lines(counts12_merge_ord[,1],line12,type="l",col="#85BF67",lwd=lwd_par)
#points(counts3_merge_ord[,1], counts3_merge_ord[,2],pch=19,col="red",lwd=lwd_par,cex=0.4)

lines(counts13_merge_ord[,1],line13,type="l",col="#D17BE9",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)

lines(counts14_merge_ord[,1],line14,type="l",col="#D17BE9",lwd=lwd_par)
#points(counts2_merge_ord[,1], counts2_merge_ord[,2],pch=19,col="blue",lwd=lwd_par,cex=0.4)

lines(counts15_merge_ord[,1],line15,type="l",col="#F15A92",lwd=lwd_par)
#points(counts3_merge_ord[,1], counts3_merge_ord[,2],pch=19,col="red",lwd=lwd_par,cex=0.4)

lines(counts16_merge_ord[,1],line16,type="l",col="#F15A92",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)

lines(counts17_merge_ord[,1],line17,type="l",col="#375F27",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)

lines(counts18_merge_ord[,1],line18,type="l",col="#375F27",lwd=lwd_par)
#points(counts4_merge_ord[,1], counts4_merge_ord[,2],pch=19,col="darkgreen",lwd=lwd_par,cex=0.4)

abline(v= center_coord + resolution/2,lty=2,col="gray")
#legend("topright",c("0818","1209","H9","WT33"), col=c("blue","red","purple","darkgreen"),inset=0.03,lwd=lwd_par,cex=0.7)

#legend("top",c("Astrocytes-B1","Astrocytes-B2","H9-D0-B1","H9-D0-B2","H9-D2-B1","H9-D2-B2","H9-D10-B1","H9-D10-B2","H9-D28-B1","H9-D28-B2","SL-D0-B1","SL-D0-B2","SL-D2-B1","SL-D2-B2","SLC-D0-B1","SLC-D0-B2","SLC-D2-B1","SLC-D2-B2"), col=c("#4E6AE1","#4E6AE1","#D5AAE1","#D5AAE1","#D17BE9","#D17BE9","#740792","#740792","#58066F","#58066F","#FF8CB8","#FF8CB8","#F15A92","#F15A92","#85BF67","#85BF67","#375F27","#375F27"),inset=0.03,lwd=lwd_par,cex=0.5)
legend("top",c("Astrocytes","H9-D0","H9-D2","H9-D10","H9-D28","SL-D0","SL-D2","SLC-D0","SLC-D2"), col=c("#4E6AE1","#D5AAE1","#D17BE9","#740792","#58066F","#FF8CB8","#F15A92","#85BF67","#375F27"),inset=0.03,lwd=lwd_par,cex=0.5)
