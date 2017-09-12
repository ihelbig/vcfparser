# vcfparser
```
setwd("~/Desktop")


library(data.table)
library(readr)
library(dplyr)

infile <- "POD-0230_POD-0230-01.trio.com.filtered.ad.de.nm.snpeff.filtered_myanno.hg19_multianno.txt"

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

table(anno_file$Other_col11)

anno_file <- subset(anno_file,anno_file$Other_col9=="PASS")
NROW(anno_file)

anno_file <- subset(anno_file,anno_file$Other_col11 == "GT:AD:DP:GQ:PL" | anno_file$Other_col11 == "GT:AD:DP:GQ:PGT:PID:PL")
NROW(anno_file)

start_time <- Sys.time()
b <- hyperloop(anno_file,"Other_col11",99,"Other_col12",100,"son")
b <- hyperloop(b,"Other_col11",99,"Other_col13",101,"mom")
b <- hyperloop(b,"Other_col11",99,"Other_col14",102,"dad")
end_time <- Sys.time()

end_time - start_time

anno_file <- b
View(anno_file[,c(1:10,80:NCOL(anno_file))])

columns_to_show <- c(1:10,100:126)

#Trio 2 son

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
               & (DP_mom >=10 & DP_son >=10 & DP_dad >=10) 
               & (GQ_mom >= 30 & GQ_son >= 30 & GQ_dad >=30)
)
NROW(rare)


# 1 - de novo
denovo_son <- subset(anno_file,
                (Other_col9 == "PASS")     
                 & (GT_son == "0/1" & GT_mom == "0/0" & GT_dad == "0/0") 
                 & (is.na(genomicSuperDups)==TRUE)
                 & (esp6500siv2_all = 0 | is.na(esp6500siv2_all)==TRUE)
                 & (anno_file[,15] = 0 | is.na(anno_file[,15])==TRUE)
                 & (ExAC_nonpsych_ALL = 0 | is.na(ExAC_nonpsych_ALL)==TRUE)
                 & (GME_AF = 0 | is.na(GME_AF)==TRUE)
                 & !(Func.refGene %in% refGene_exclude)
                 & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                 
                 & (DP_mom >=10 & DP_son >=10 & DP_dad >=10)
                 & (GQ_mom >= 30 & GQ_son >=30 & GQ_dad >=30)
)
NROW(denovo_son)
View(denovo_son[,columns_to_show])

# 2 - recessive 
hom_son <- subset(anno_file,
                  (Other_col9 == "PASS")     
                  & (GT_son == "1/1" & GT_mom == "0/1" & GT_dad == "0/1") 
                  & (is.na(genomicSuperDups)==TRUE)
                  & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
                  & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
                  & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
                  & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
                  & !(Func.refGene %in% refGene_exclude)
                  & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                  
                  & (DP_mom >=10 & DP_son >=10 & DP_dad >=10)
                  & (GQ_mom >= 30 & GQ_son >=30 & GQ_dad >=30)
)
NROW(hom_son)
View(hom_son[,columns_to_show])

# 3 - compound
# 3 - compound
mat_son <- subset(anno_file,
                       (Other_col9 == "PASS")     
                       & (GT_son == "0/1" & GT_mom == "0/1" & GT_dad == "0/0") 
                       & (is.na(genomicSuperDups)==TRUE)
                       & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
                       & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
                       & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
                       & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
                       & !(Func.refGene %in% refGene_exclude)
                       & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                       
                       & (DP_mom >=10 & DP_son >=10 & DP_dad >=10)
                       & (GQ_mom >= 30 & GQ_son >=30 & GQ_dad >=30)
)
NROW(mat_son)

pat_son <- subset(anno_file,
                       (Other_col9 == "PASS")     
                       & (GT_son == "0/1" & GT_mom == "0/0" & GT_dad == "0/1") 
                       & (is.na(genomicSuperDups)==TRUE)
                       & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
                       & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
                       & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
                       & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
                       & !(Func.refGene %in% refGene_exclude)
                       & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                       
                       & (DP_mom >=10 & DP_son >=10 & DP_dad >=10)
                       & (GQ_mom >= 30 & GQ_son >=30 & GQ_dad >=30)
)
NROW(pat_son)

#combine
matgene_son <- data.frame(table(mat_son$Gene.refGene)); names(matgene_son) <- c("gene","mat")
patgene_son <- data.frame(table(pat_son$Gene.refGene)); names(patgene_son) <- c("gene","pat")
compound_son <- merge(matgene_son,patgene_son, by.x="gene",by.y="gene",all.x=T, all.y=T)

compound2_son <- subset(compound_son,
                             (compound_son$mat < 3 & is.na(compound_son$mat)==FALSE) 
                             & (compound_son$pat < 3 & is.na(compound_son$mat)==FALSE)
)
NROW(compound2_son)
comp_list_son <- as.vector(compound2_son$gene)

compound_het_son <- subset(rare,
                                (GT_son == "0/1")
                                & (Gene.refGene %in% comp_list_son)
                                & !(Gene.refGene %in% Genes_exclude)
)

View(compound_het_son[,columns_to_show])

# X - linked
X_son <- subset(anno_file,
                (Other_col9 == "PASS")     
                & (GT_son == "1/1" & GT_mom == "0/1" & GT_dad == "0/0") 
                & (Chr=="X")
                & (is.na(genomicSuperDups)==TRUE)
                & (esp6500siv2_all <= 0.01 | is.na(esp6500siv2_all)==TRUE)
                & (anno_file[,15] <= 0.01 | is.na(anno_file[,15])==TRUE)
                & (ExAC_nonpsych_ALL <= 0.01 | is.na(ExAC_nonpsych_ALL)==TRUE)
                & (GME_AF <= 0.01 | is.na(GME_AF)==TRUE)
                & !(Func.refGene %in% refGene_exclude)
                & !(ExonicFunc.refGene %in% ExonicFunc_exclude)
                
                & (DP_mom >=10 & DP_son >=10 & DP_dad >=10)
                & (GQ_mom >= 30 & GQ_son >=30 & GQ_dad >=30)
)
NROW(X_son)
View(X_son[,columns_to_show])



```
