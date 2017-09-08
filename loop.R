```
a <- a[,1:81]
colx <- NCOL(a)
a[,(colx+1):(colx+21)] <- NA
namo <- names(a)
namo[(colx+1):(colx+21)] <- c(
  
  "GT1_GT","GT1_GQ","GT1_AD1","GT1_AD2","GT1_PL1","GT1_PL2","GT1_PL3",
  "GT2_GT","GT2_GQ","GT2_AD1","GT2_AD2","GT2_PL1","GT2_PL2","GT2_PL3",
  "GT3_GT","GT3_GQ","GT3_AD1","GT3_AD2","GT3_PL1","GT3_PL2","GT3_PL3"
)
names(a) <- namo
names(a)

start_time <- Sys.time()
for (i in 1:2000)
{
  gt1 <- gt_extract(a$Otherinfo_8[i],a$Otherinfo_9[i])
  gt2 <- gt_extract(a$Otherinfo_8[i],a$Otherinfo_10[i])
  gt3 <- gt_extract(a$Otherinfo_8[i],a$Otherinfo_11[i])
  a[i,(colx+1):(colx+7)] <- gt1
  a[i,(colx+8):(colx+14)] <- gt2
  a[i,(colx+15):(colx+21)] <- gt3
}
end_time <- Sys.time()
View(a[,70:102])
```
