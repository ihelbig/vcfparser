
```
library(tidyr)
a1 <- a[,1:81]
start_time <- Sys.time()
uniq <- unique(a1$Otherinfo_8)

ix = 1
uniq[ix]
luniq <- length(uniq)

#Generate number of dataframes that can be used

for (i in 1:luniq)
{
  bdy <- subset(a1,a1$Otherinfo_8==uniq[i])
  attach(bdy)
  into = c(strsplit(uniq[i],":")[[1]]) #check VCF first, this can differ
  ix = 1
  into1 <- paste(into,ix,sep='_')
  ad1 <- paste(c("AD1","AD2"),ix,sep='_')
  pl1 <- paste(c("PL1","PL2","PL3"),ix,sep='_')
  bdy <- separate(bdy,Otherinfo_9,into1, sep =":",remove=F) #tidyr function
  bdy <- separate(bdy,AD_1,ad1, sep =",",remove=F) #tidyr function
  bdy <- separate(bdy,PL_1,pl1, sep =",",remove=F) #tidyr function
  
  ix = 2
  into2 <- paste(into,ix,sep='_')
  ad2 <- paste(c("AD1","AD2"),ix,sep='_')
  pl2 <- paste(c("PL1","PL2","PL3"),ix,sep='_')
  bdy <- separate(bdy,Otherinfo_10,into2, sep =":",remove=F) #tidyr function
  bdy <- separate(bdy,AD_2,ad2, sep =",",remove=F) #tidyr function
  bdy <- separate(bdy,PL_2,pl2, sep =",",remove=F) #tidyr function
  
  ix = 3
  into3 <- paste(into,ix,sep='_')
  ad3 <- paste(c("AD1","AD2"),ix,sep='_')
  pl3 <- paste(c("PL1","PL2","PL3"),ix,sep='_')
  bdy <- separate(bdy,Otherinfo_11,into3, sep =":",remove=F) #tidyr function
  bdy <- separate(bdy,AD_3,ad3, sep =",",remove=F) #tidyr function
  bdy <- separate(bdy,PL_3,pl3, sep =",",remove=F) #tidyr function
  
  assign(paste("gt_array",i,sep=""),bdy)
}

for(i in 1:luniq)
{
  x <- get(ls(pattern="^gt_array")[i])
  print(as.numeric(x$DP_1[1:5]))
}
names(gt_array1)

View(bdy[,79:109])

xlen = length(ls(pattern="^gt_array"))

# create ordered array a[,c(1:78,"Otherinfo_9","Otherinfo10"...]
namex <- names(a1[,1:78])
orderx <- c(
  "Otherinfo_9","Otherinfo_10","Otherinfo_11",
  "GT_1", "GT_2", "GT_3",
  "AD1_1", "AD2_1", "AD1_2", "AD2_2", "AD1_3", "AD2_3",
  "DP_1", "DP_2", "DP_3",
  "GQ_1", "GQ_2","GQ_3",
  "PL1_1", "PL2_1", "PL3_1", 
  "PL1_2", "PL2_2", "PL3_2",
  "PL1_3", "PL2_3", "PL3_3")  
combinedx <- c(namex,orderx)
colxx <- length(combinedx)
df_final <- as.data.frame((matrix(ncol=colxx,nrow=0)))
names(df_final) <- combinedx

for(i in 1:xlen)
{

dff <- get(ls(pattern="^gt_array")[i])
dff1 <- dff[,1:78]
dff2 <- dff[,orderx]
dffx <- cbind(dff1,dff2) # this creates a dataframe in the same order for every column
df_final <- rbind(df_final,dffx)
}

df_final <- df_final[order( df_final$Chr, df_final$Start ),]
rownames(df_final) <- NULL
end_time <- Sys.time()
end_time - start_time
View(df_final[,1:10,70:105])

#Time difference of 2.498291 secs instead of 34 min during the regular for... loop
```
