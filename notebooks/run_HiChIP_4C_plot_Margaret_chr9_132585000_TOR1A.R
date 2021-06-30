setwd("/Users/lauradonohue/Khavari_Bioinformatics/Project/Wernig_collab/HiChIP/4C_plots_021720")
#setwd("/Users/mguo123/Documents/pan_omics_psych/notebooks")
library(Hmisc)
library(zoo)
library(tidyverse)

######################################
# Setup
######################################
# read in dataframes and set save path (FILL IN Accordingly)
save_dir = "/Users/lauradonohue/Khavari_Bioinformatics/Project/Wernig_collab/HiChIP/4C_plots_021720/"
dir_4C = "/Users/lauradonohue/Khavari_Bioinformatics/Project/Wernig_collab/HiChIP/4C_plots_021720/"
color_df  = read.csv("colors.csv", row.names=1)
color_arr = color_df$color_id
names(color_arr) = rownames(color_df)
valid_stat_df = read.table(paste0(dir_4C,"valid_stat.tsv"),sep='\t',
                           row.names=1,header=T)

# color_df  = read.csv("colors.csv", stringsAsFactor=F,row.names=1)
# color_arr = color_df$color_id
# names(color_arr) = rownames(color_df)

# set variables (FILL IN PER RUN)
samples = c('Astro_B1', 'Astro_B2', 'H9-B1','H9-B2','H9D2-B1','H9D2-B2','H9D10-B1','H9D10-B2',
            'H9D28-B1','H9D28-B2','SL-B1','SL-B2','SLC-B1','SLC-B2','SLD2-B1','SLD2-B2','SLCD2-B1','SLCD2-B2')
gene_sym="TOR1A"
chrom=9
position=132585000 # on file
center_coord=132585000
resolution=5000

# derived variables
### DEBUG if the file reading in doesn't work by 
# list.files(dir_4C,pattern='*bed')
file_prefix = paste0(dir_4C,'4C_chr',chrom, '_', position, '_', gene_sym, '_')
file_suffix = paste0('_',as.character(resolution), '.matrix.bed')
save_file = paste0(save_dir,gene_sym, '_4tissues.pdf')
print(save_file)


######################################
# Reading in files
######################################

## Define reading table function
# sample : <str> REQUIRED
# norm_method: <str> 'valid', 'sum', default sets table_norm = 1
# rollmean_win: <int> Pick your smoothing window (1 means no smoothing, 2 averages over 2 adjacent bins, etc. )
read_bed = function(sample, norm_method='valid', rollmean_win =1){
  # 1. read bed
  tissue = strsplit(sample, '[_-]')[[1]][1]
  table = read.table(paste0(file_prefix, sample, file_suffix))
  # 2. normalize
  if (norm_method=='valid'){ # from # valid interactions
    table_norm = valid_stat_df[sample, "X"]
  }
  else if (norm_method=='sum'){ # sum/1000
    table_norm = sum(table[,7])/1e3
  }
  else{ # no norm
    table_norm = 1
  }
  
  # Process table
  counts_left=table[table[,2]== center_coord,c(5,7)]
  counts_right=table[table[,5]== center_coord,c(2,7)]
  colnames(counts_left) = c(1,2)
  colnames(counts_right) = c(1,2)
  counts_merge=rbind(counts_left, counts_right)
  counts_merge_ord = counts_merge[order(counts_merge[,1]),]
  counts_merge_ord[,2] = counts_merge_ord[,2]/table_norm
  counts_merge_ord[,1] = counts_merge_ord[,1] + resolution/2
  line = rollmean(counts_merge_ord[,2], rollmean_win,fill=c(0,0,0))
  
  
  df = data.frame(x = counts_merge_ord[,1], y= line, sample = sample, tissue = tissue)
  
  return (df)
}

# read in chosen samples
tables = lapply(samples, read_bed)
table_all = bind_rows(tables)
table_all_avg = table_all%>%
  group_by(x,tissue)%>%
  summarise(y=mean(y))


######################################
# Plot
######################################
# plotting variables
xlim_minus=7e4 # bp upstream of anchor
xlim_plus=1e5 # bp downstream of anchor
xlim_par= c((center_coord)-xlim_minus,(center_coord)+xlim_plus)
ylim_par=c(-0.01,0.2) # y axis limits
lwd_par = .5 # line thickness #### NOT USED
xlab="Position on chr3"
ylab="Filtered reads"
plot_name=paste("Virtual 4C anchor 5kb bin for gene: ",gene_sym, ' at chr',chrom, ':',center_coord,sep="")


# actual plot 
ggplot(table_all_avg, aes(x=x, y=y, color=tissue))+geom_line(lwd=lwd_par)+  
  scale_colour_manual(values=c(Astro="#4E6AE1",H9="#D5AAE1",H9D2="#D17BE9",H9D10="#740792",H9D28="#58066F",SL="#FF8CB8",SLD2="#F15A92",SLC="#85BF67",SLCD2="#375F27"))+
  labs(title=plot_name,
        x =xlab, y = ylab)+
  xlim(xlim_par) + ylim(ylim_par)+
  theme_classic()+ theme(plot.title = element_text(size=10))+
  geom_vline(xintercept = center_coord + resolution/2,color='gray',lty=2)+
  geom_vline(xintercept = 132565943,color='red',lty=2)+
  geom_vline(xintercept = 132596899,color='red',lty=2)+
  geom_vline(xintercept = 132597721,color='red',lty=2)+
  geom_vline(xintercept = 132597840,color='red',lty=2)+
  geom_vline(xintercept = 132599101,color='red',lty=2)+
  geom_vline(xintercept = 132601229,color='red',lty=2)

# saving plot
ggsave(save_file) # CAN CHANGE width and height
