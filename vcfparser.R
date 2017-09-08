gt_extract <- function(format,gt_field)
{
  colx <- strsplit(as.character(format),":")[[1]]
  gtx <- as.data.frame(matrix(ncol=length(colx),nrow=1))
  colnames(gtx) <- colx
  gt1 <- as.character(gt_field)
  if (!is.na(gt1)) {gtx[1,] <- strsplit(gt1,":")[[1]]}

  gtx2 <- gtx
  num <- ncol(gtx2)
  gtx2[,num+1] <- NA #AD1
  gtx2[,num+2] <- NA #AD2
  gtx2[,num+3] <- NA #PL1
  gtx2[,num+4] <- NA #PL2
  gtx2[,num+5] <- NA #PL3
  
  namo <- names(gtx2)
  namo[num+1] <- "AD1"
  namo[num+2] <- "AD2"
  namo[num+3] <- "PL1"
  namo[num+4] <- "PL2"
  namo[num+5] <- "PL3"
  colnames(gtx2) <- namo
  
  #First line 
  #line any line (linex)
  
  ad <- as.character(gtx2$AD[1])
  pl <- as.character(gtx2$PL[1])
  if (!is.na(ad)) {ad_s <- strsplit(ad,",")[[1]]; gtx2[1,(num+1):(num+2)] <- ad_s}
  if (!is.na(pl)) {pl_s <- strsplit(pl,",")[[1]]; gtx2[1,(num+3):(num+5)] <- pl_s}
  
  table_gt <- gtx2[,c("GT","GQ","AD1","AD2","PL1","PL2","PL3")]
  
  return(table_gt)
  # table_gt$GT
  # table_gt$GQ
  # table_gt$AD1
  # table_gt$AD2
  # table_gt$PL1
  # table_gt$PL2
  # table_gt$PL3
}
