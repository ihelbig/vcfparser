# vcfparser
```
setwd("~/Desktop")


library(data.table)
library(readr)
library(dplyr)

infile <- "POD-0241_POD-0241-01.trio.com.filtered.ad.de.nm.snpeff.filtered_myanno.hg19_multianno.txt"

file_header <- read.table(infile,nrows=1)
file_body <- read.table(infile,nrows=1,skip=1,sep="\t")
file_header_l <- length(file_header)
file_body_l <- length(file_body)

file_header_l
file_body_l

# diff_l is the difference in colums between header and file
#130 columns in combined file for trio

diff_l <- file_body_l - file_header_l

new_header1 <- as.vector(as.matrix(file_header[1,]))
new_header2 <- paste("Other_col",rep(1:diff_l),sep="")
new_header3 <- c(new_header1,new_header2)

#Manually confirming that all column header start with symbol, replacing 1000G headers
#new_header3[14] <- "X1000g2015aug_all"
#new_header3[15] <- "X1000g2015aug_eur"


anno_file <- read_tsv(infile,skip=1,na=".")
colnames(anno_file) <- new_header3

#Precleaning

anno_file <- subset(anno_file,anno_file$Other_col9=="PASS")
NROW(anno_file)

anno_file <- subset(anno_file,anno_file$Other_col11 == "GT:AD:DP:GQ:PL" | anno_file$Other_col11 == "GT:AD:DP:GQ:PGT:PID:PL")
NROW(anno_file)

start_time <- Sys.time()
b <- hyperloop(anno_file,"Other_col11",99,"Other_col12",100,"daughter")
b <- hyperloop(b,"Other_col11",99,"Other_col13",101,"mom")
b <- hyperloop(b,"Other_col11",99,"Other_col14",102,"dad")
end_time <- Sys.time()

end_time - start_time

anno_file <- b
View(anno_file[,c(1:10,80:NCOL(anno_file))])

columns_to_show <- c(1:10,100:126)

#rare variants
rare <- subset(anno_file,
               (Other_col9 == "PASS")
               & (is.na(genomicSuperDups)==TRUE)
               & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
               & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
               & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
               & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
               & !(Func.refGene %in% refGene_exclude)
               & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
               & (DP_mom >=10 & DP_daughter >=10 & DP_dad >=10) 
               & (GQ_mom >= 30 & GQ_daughter >= 30 & GQ_dad >=30)
)
NROW(rare)


# 1 - de novo
denovo_daughter <- subset(anno_file,
                          (Other_col9 == "PASS")     
                          & (GT_daughter == "0/1" & GT_mom == "0/0" & GT_dad == "0/0") 
                          & (is.na(genomicSuperDups)==TRUE)
                          & (esp6500siv2_all = 0 | is.na(esp6500siv2_all)==TRUE)
                          & (anno_file[,15] = 0 | is.na(anno_file[,15])==TRUE)
                          & (ExAC_nonpsych_ALL = 0 | is.na(ExAC_nonpsych_ALL)==TRUE)
                          & (GME_AF = 0 | is.na(GME_AF)==TRUE)
                          & !(Func.refGene %in% refGene_exclude)
                          & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                          
                          & (DP_mom >=10 & DP_daughter >=10 & DP_dad >=10)
                          & (GQ_mom >= 30 & GQ_daughter >=30 & GQ_dad >=30)
)
NROW(denovo_daughter)
View(denovo_daughter[,columns_to_show])

# 2 - recessive 
hom_daughter <- subset(anno_file,
                       (Other_col9 == "PASS")     
                       & (GT_daughter == "1/1" & GT_mom == "0/1" & GT_dad == "0/1") 
                       & (is.na(genomicSuperDups)==TRUE)
                       & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
                       & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
                       & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
                       & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
                       & !(Func.refGene %in% refGene_exclude)
                       & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                       
                       & (DP_mom >=10 & DP_daughter >=10 & DP_dad >=10)
                       & (GQ_mom >= 30 & GQ_daughter >=30 & GQ_dad >=30)
)
NROW(hom_daughter)
View(hom_daughter[,columns_to_show])

# 3 - compound
# 3 - compound
mat_daughter <- subset(anno_file,
                       (Other_col9 == "PASS")     
                       & (GT_daughter == "0/1" & GT_mom == "0/1" & GT_dad == "0/0") 
                       & (is.na(genomicSuperDups)==TRUE)
                       & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
                       & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
                       & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
                       & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
                       & !(Func.refGene %in% refGene_exclude)
                       & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                       
                       & (DP_mom >=10 & DP_daughter >=10 & DP_dad >=10)
                       & (GQ_mom >= 30 & GQ_daughter >=30 & GQ_dad >=30)
)
NROW(mat_daughter)

pat_daughter <- subset(anno_file,
                       (Other_col9 == "PASS")     
                       & (GT_daughter == "0/1" & GT_mom == "0/0" & GT_dad == "0/1") 
                       & (is.na(genomicSuperDups)==TRUE)
                       & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
                       & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
                       & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
                       & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
                       & !(Func.refGene %in% refGene_exclude)
                       & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                       
                       & (DP_mom >=10 & DP_daughter >=10 & DP_dad >=10)
                       & (GQ_mom >= 30 & GQ_daughter >=30 & GQ_dad >=30)
)
NROW(pat_daughter)

#combine
matgene_daughter <- data.frame(table(mat_daughter$Gene.refGene)); names(matgene_daughter) <- c("gene","mat")
patgene_daughter <- data.frame(table(pat_daughter$Gene.refGene)); names(patgene_daughter) <- c("gene","pat")
compound_daughter <- merge(matgene_daughter,patgene_daughter, by.x="gene",by.y="gene",all.x=T, all.y=T)

compound2_daughter <- subset(compound_daughter,
                             (compound_daughter$mat < 3 & is.na(compound_daughter$mat)==FALSE) 
                             & (compound_daughter$pat < 3 & is.na(compound_daughter$mat)==FALSE)
)
NROW(compound2_daughter)
comp_list_daughter <- as.vector(compound2_daughter$gene)

compound_het_daughter <- subset(rare,
                                (GT_daughter == "0/1")
                                & (Gene.refGene %in% comp_list_daughter)
                                & !(Gene.refGene %in% Genes_exclude)
)
View(compound_het_daughter[,columns_to_show])

```
