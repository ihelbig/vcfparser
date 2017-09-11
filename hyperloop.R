
```
hyperloop <- function(vcf, format, pos_format, gt_field, pos_gt_field, suffix)
  
{
  library(tidyr)
  
  start_time <- Sys.time()
  uniq <- unique(vcf[[format]])
  print(uniq)
  luniq <- length(uniq)
  
  #Generate number of dataframes that can be used
  
  for (i in 1:luniq)
  {
    bdy <- subset(vcf,vcf[[format]]==uniq[i])
    col_begin <- ncol(bdy)
    attach(bdy)
    into = c(strsplit(uniq[i],":")[[1]]) #check VCF first, this can differ
    into <- paste(into,suffix,sep="_")
    ad1 <- paste(c("AD1","AD2"),suffix,sep="_")
    pl1 <- paste(c("PL1","PL2","PL3"),suffix,sep="_")
    bdy <- separate_(bdy,gt_field,into, sep =":",remove=F) #tidyr function
    bdy <- separate_(bdy,paste("AD",suffix,sep="_"),ad1, sep =",",remove=F) #tidyr function
    bdy <- separate_(bdy,paste("PL",suffix,sep="_"),pl1,sep =",",remove=F) #tidyr function
    col_end <- ncol(bdy)
    new_cols <- col_end - col_begin
    additional_gt <- col_begin - pos_gt_field
    orderx <- c(1:pos_gt_field,((col_end-additional_gt)+1):col_end,(pos_gt_field+1):(col_end-additional_gt))
    bdy <- bdy[,orderx]
    names(bdy)
    assign(paste("gt_array",i,sep=""),bdy)
  }

  xlen = length(ls(pattern="^gt_array"))
  
  #
  namex <- names(vcf)
  #orderx <- c("GT", "AD1", "AD2", "DP", "GQ", "PL1", "PL2", "PL3")
  orderx <- paste(c("GT", "AD1", "AD2", "DP", "GQ", "PL1", "PL2", "PL3"),"_",suffix,sep='')
  combinedx <- c(namex,orderx)
  colxx <- length(combinedx)
  df_final <- as.data.frame((matrix(ncol=colxx,nrow=0)))
  names(df_final) <- combinedx
  print(combinedx)
  
  for(i in 1:xlen)
  {
    dff <- get(ls(pattern="^gt_array")[i])

    dff1 <- dff[,1:col_begin]
  
    dff2 <- dff[,orderx]
   
    dffx <- cbind(dff1,dff2) # this creates a dataframe in the same order for every column
   
    df_final <- rbind(df_final,dffx)
 
  }
  
  df_final <- df_final[order( df_final$Chr, df_final$Start ),]
  rownames(df_final) <- NULL
  end_time <- Sys.time()
  print(end_time - start_time)

  return(df_final)

}
```
