
tab = read.csv(textConnection("
GPU, MLUPS, MLUPS_pf
V100  -  Gadi, 1599.3, 641.4
A100  -  Athena, 2857.6, 1176.9
A100 DGX  -  Gadi, 3450.2, 1353.0
H100  -  Bunya, 3737.8, 1526.0
MI250  -  Setonix, 1754.8, 352.0
"),stringsAsFactor=FALSE)

pdf("perf_cumulant.pdf",width=8, height = 5)
par(mfrow=c(1,2),mar=c(0.4,4.1,4.1,1.1))
ret = barplot(tab$MLUPS,ylab="MLUps",main="d3q27_cumulant")
text(ret,tab$MLUPS/2,labels = tab$GPU,srt=90)
ret = barplot(tab$MLUPS_pf,ylab="MLUps",main="d3q27_pf_velocity")
text(ret,tab$MLUPS_pf/2,labels = tab$GPU,srt=90)
dev.off()
