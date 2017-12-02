setwd('~/Desktop/new couses/Genetic Analysis/Assignment/final project/result data')

load('Yu Wang_11.Rdata')
load('Yu Wang_12.Rdata')
load('Yu Wang_13.Rdata')
load('Yu Wang_21.Rdata')
load('Yu Wang_22.Rdata')
load('Yu Wang_23.Rdata')
load('D10h1er0_05.RData')
load('results.matrix.RData')
load('Huiyin_D10h1M0.05.RData')
load('Huiyin_D10h1M0.2.RData')
Coverage.1<- c(0.248, 0.414, 0.855, 0.896, 0.777, 0.449, 0.407, 0.452, 0.718, 0.818, 0.729, 0.460, 0.437, 0.463,
               0.738, 0.819, 0.744, 0.444, 0.415, 0.428, 0.743, 0.830, 0.739, 0.460, 0.451, 0.482, 0.794, 0.910,
               0.847, 0.386, 0.265)
Coverage.2<- c(0.116, 0.230, 0.744, 0.977, 0.958, 0.949, 0.952, 0.928, 0.928, 0.982, 0.970, 0.972, 0.960, 0.947,
               0.947, 0.984, 0.945, 0.940, 0.949, 0.942, 0.959, 0.983, 0.931, 0.936, 0.953, 0.948, 0.963, 0.975,
               0.715, 0.231, 0.122)
Size.1<-c(6.843,  7.095,  7.970,  8.509,  8.873,  9.041,  8.934,  9.741, 10.304, 10.596, 10.187, 10.052,
          9.780, 10.194, 10.188, 10.299, 10.428, 10.005,  9.617,  9.897, 10.499, 10.597, 10.021, 9.788,
          8.916 , 9.045 , 8.770,  8.250 , 7.707 , 7.268 , 6.638)
Size.2<-c(16.49675, 17.92035, 18.81162, 19.70362,20.02545, 20.09832, 20.15865, 22.28095, 23.46287,
                 24.37175, 24.50705, 24.03730, 22.84010, 24.59960, 25.18600, 25.64740, 25.27442, 23.87430,
                 22.54930, 23.84570, 24.34485, 24.12022, 23.17627, 21.74877, 19.40365, 19.62382 ,19.60897,
                 19.25750, 18.46410, 17.35625, 16.07552)
Distence<-c(4.8490, 3.8406, 3.1670, 2.8220, 4.1006, 4.7514, 4.8230 ,5.1516, 4.9016, 4.5770, 5.0122, 5.3442,
            5.1620, 5.2332, 4.8248, 4.4000, 4.6790, 5.1390, 5.2000, 5.2690, 4.8168, 4.3230, 4.7878, 5.0990,
            4.7980, 4.4830, 3.7220, 2.7710, 2.9264, 3.6912, 4.8780)
result.p13er0.2<-cbind(Coverage.1,Coverage.2,Size.1,Size.2,Distence)
result.p13M0.05<-cbind(Coverage.1,Coverage.2,Size.1,Size.2,Distence)
result.p13M0.2<-cbind(Coverage.1,Coverage.2,Size.1,Size.2,Distence)
rm(CI.size1,CI.size2,CIs.1,CIs.2,Coverage.1,Coverage.2,Distence,Ef.size,Effect,i,         
   Indicator.1,Indicator.2,j,lb,map,Map.dens,n,N,
   Posi, Q.dist1,Chrom.length,sim.data,Size.1,Size.2,temp1,tt1,ub)
Posi <- c(0,1.7,3.3,5,6.7,8.3,10,11.7,13.3,15,16.7,18.3,20,21.7,23.3,
          25,26.7,28.3,30,31.7,33.3,35,36.7,38.3,40,41.7,43.3,45,46.7,48.3,
          50)
Ef_p_ld10<-cbind(result.p11[,1],result.p12[,1],result.p13[,1])
Ef_p_bd10<-cbind(result.p11[,2],result.p12[,2],result.p13[,2])
Ef_s_ld10<-cbind(result.p11[,3],result.p12[,3],result.p13[,3])
Ef_s_bd10<-cbind(result.p11[,4],result.p12[,4],result.p13[,4])

Ef_p_ld5<-cbind(result.p21[,1],result.p22[,1],result.p23[,1])
Ef_p_bd5<-cbind(result.p21[,2],result.p22[,2],result.p23[,2])
Ef_s_ld5<-cbind(result.p21[,3],result.p22[,3],result.p23[,3])
Ef_s_bd5<-cbind(result.p21[,4],result.p22[,4],result.p23[,4])

Dis_D10h1_er<-cbind(result.p13[,5],result.p13er0.05[,5],result.p13er0.2[,5])

Dis_D10h1_mr<-cbind(result.p13[,5],result.p13M0.05[,5],result.p13M0.2[,5])

##############################################################
# Map.dens <- c(6,11)                                        #
# Ef.size <-c(0.05,0.5,1)                                    #
# p11：Marker number 6(10cM/each), effect size 0.05          #
# p21: Marker number 11(5cM/each), effect size 0.05          #
# p12：Marker number 6(10cM/each), effect size 0.5           #
# p22: Marker number 11(5cM/each), effect size 0.5           #
# p13：Marker number 6(10cM/each), effect size 1             #
# p23: Marker number 11(5cM/each), effect size 1             #
# 1: LOD method, 2: Boostrap method                          #
##############################################################

matplot(Posi,result.p11[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        main='Delta=10,h=0.05',font=2,font.lab=2)
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0.8,0.8,0.8,0.8,0.8,0.8),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

matplot(Posi,result.p11[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(47,50),font=2,font.lab=2,
        main='Delta=10,h=0.05')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(45,45,45,45,45,45),c(0,10,20,30,40,50),c(49.4,49.4,49.4,49.4,49.4,49.4),
         lty = 4,col=6)

matplot(Posi,result.p21[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Delta=5,h=0.05')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8),
         c(0,5,10,15,20,25,30,35,40,45,50),c(1,1,1,1,1,1,1,1,1,1,1),lty = 4,col=6)

matplot(Posi,result.p21[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(47,50),font=2,font.lab=2,
        main='Delta=5,h=0.05')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(45,45,45,45,45,45,45,45,45,45,45),
         c(0,5,10,15,20,25,30,35,40,45,50),c(49,49,49,49,49,49,49,49,49,49,49),lty = 4,col=6)

matplot(Posi,result.p13[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Delta=10,h=1')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0.8,0.8,0.8,0.8,0.8,0.8),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

matplot(Posi,result.p13[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(0,14),font=2,font.lab=2,
        main='Delta=10,h=1')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(11,11,11,11,11,11),
         lty = 4,col=6)

matplot(Posi,result.p23[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Delta=5,h=1')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8),
         c(0,5,10,15,20,25,30,35,40,45,50),c(1,1,1,1,1,1,1,1,1,1,1),lty = 4,col=6)


matplot(Posi,result.p23[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(0,14),font=2,font.lab=2,
        main='Delta=5,h=1')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(0,0,0,0,0,0,0,0,0,0,0),
         c(0,5,10,15,20,25,30,35,40,45,50),c(8,8,8,8,8,8,8,8,8,8,8),lty = 4,col=6)


matplot(Posi,result.p12[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Delta=10,h=0.5')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0.8,0.8,0.8,0.8,0.8,0.8),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

matplot(Posi,result.p12[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(5,35),font=2,font.lab=2,
        main='Delta=10,h=0.5')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(5,5,5,5,5,5),c(0,10,20,30,40,50),c(25,25,25,25,25,25),
         lty = 4,col=6)

matplot(Posi,result.p22[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Delta=5,h=0.5')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8),
         c(0,5,10,15,20,25,30,35,40,45,50),c(1,1,1,1,1,1,1,1,1,1,1),lty = 4,col=6)


matplot(Posi,result.p22[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(5,35),font=2,font.lab=2,
        main='Delta=5,h=0.5')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(5,5,5,5,5,5,5,5,5,5,5),
         c(0,5,10,15,20,25,30,35,40,45,50),c(25,25,25,25,25,25,25,25,25,25,25),
         lty = 4,col=6)
##################################### Er

matplot(Posi,result.p13er0.05[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.5,1.15),
        font=2,font.lab=2,main='Delta=10,h=1,error rate=0.05')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0.5,0.5,0.5,0.5,0.5,0.5),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

matplot(Posi,result.p13er0.05[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(0,20),font=2,font.lab=2,
        main='Delta=10,h=1,error rate=0.05')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(11,11,11,11,11,11),
         lty = 4,col=6)

matplot(Posi,result.p13er0.2[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0,1.4),
        font=2,font.lab=2,main='Delta=10,h=1,error rate=0.2')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

matplot(Posi,result.p13er0.2[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(0,35),font=2,font.lab=2,
        main='Delta=10,h=1,error rate=0.2')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(22,22,22,22,22,22),
         lty = 4,col=6)

########################### Missing marker genotype analysis.###########

matplot(Posi,result.p13M0.05[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Delta=10,h=1,Missing rate=0.05')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0.85,0.85,0.85,0.85,0.85,0.85),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

matplot(Posi,result.p13M0.05[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(0,14),font=2,font.lab=2,
        main='Delta=10,h=1,Missing rate=0.05')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(8,8,8,8,8,8),
         lty = 4,col=6)

## MR=0.2
matplot(Posi,result.p13M0.2[,1:2],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Delta=10,h=1,Missing rate=0.2')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0.85,0.85,0.85,0.85,0.85,0.85),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

matplot(Posi,result.p13M0.2[,3:4],type = 'o',col = c('black','blue'),pch = c(6,8),
        ylab = 'CI size(cM)',xlab = 'QTL(cM)',ylim=c(0,14.5),font=2,font.lab=2,
        main='Delta=10,h=1,Missing rate=0.2')
legend('topright',legend = c('1 LOD method','Boostrap method'),col = c('black','blue'),
       lty = 1:2,pch = c(6,8))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(8,8,8,8,8,8),
         lty = 4,col=6)

########## Directly comparision ################

matplot(Posi,Ef_p_ld10,type = 'o',col = c('black','blue','red'),pch = c(6,8,15),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Lod method, Delta=10')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('h=0.05','h=0.5','h=1'),col = c('black','blue','red'),
       lty = 1:2,pch = c(6,8,15))
segments(c(0,10,20,30,40,50),c(0.85,0.85,0.85,0.85,0.85,0.85),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

# Lod method,Delta=5

matplot(Posi,Ef_p_ld5,type = 'o',col = c('black','blue','red'),pch = c(6,8,15),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Lod method, Delta=5')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('h=0.05','h=0.5','h=1'),col = c('black','blue','red'),
       lty = 1:2,pch = c(6,8,15))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85),
         c(0,5,10,15,20,25,30,35,40,45,50),c(1,1,1,1,1,1,1,1,1,1,1),lty = 4,col=6)


# Boostrap method, Delta=10

matplot(Posi,Ef_p_bd10,type = 'o',col = c('black','blue','red'),pch = c(6,8,15),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Boostrap method, Delta=10')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('h=0.05','h=0.5','h=1'),col = c('black','blue','red'),
       lty = 1:2,pch = c(6,8,15))
segments(c(0,10,20,30,40,50),c(0.85,0.85,0.85,0.85,0.85,0.85),c(0,10,20,30,40,50),c(1,1,1,1,1,1),
         lty = 4,col=6)

# Boostrap method, Delta=5

matplot(Posi,Ef_p_bd5,type = 'o',col = c('black','blue','red'),pch = c(6,8,15),
        ylab = 'Proportion of coverage',xlab = 'QTL(cM)',ylim = c(0.85,1.06),
        font=2,font.lab=2,main='Boostrap method, Delta=5')
abline(h=0.95,col='red',lty=3)
legend('topright',legend = c('h=0.05','h=0.5','h=1'),col = c('black','blue','red'),
       lty = 1:2,pch = c(6,8,15))
segments(c(0,5,10,15,20,25,30,35,40,45,50),c(0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85),
         c(0,5,10,15,20,25,30,35,40,45,50),c(1,1,1,1,1,1,1,1,1,1,1),lty = 4,col=6)

# Estimation Accuracy analysis##
# Error rate
matplot(Posi,Dis_D10h1_er,type = 'o',col = c('black','blue','red'),pch = c(6,8,15),
        ylab = 'Dis from TureQTL(cM)',xlab = 'QTL(cM)',ylim=c(0,8),font=2,font.lab=2,
        main='QTL Accuracy with Genotype Errer(Delt=10,h=1)')
legend('topright',legend = c('ER=0','ER=0.05','ER=0.2'),col = c('black','blue','red'),
       lty = 1:2,pch = c(6,8,15))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(5,5,5,5,5,5),
         lty = 4,col=6)


# Missing value.
matplot(Posi,Dis_D10h1_mr,type = 'o',col = c('black','blue','red'),pch = c(6,8,15),
        ylab = 'Dis from TureQTL(cM)',xlab = 'QTL(cM)',ylim=c(0,3),font=2,font.lab=2,
        main='QTL Accuracy with Genotype Missing(Delt=10,h=1)')
legend('topright',legend = c('MR=0','MR=0.05','MR=0.2'),col = c('black','blue','red'),
       lty = 1:2,pch = c(6,8,15))
segments(c(0,10,20,30,40,50),c(0,0,0,0,0,0),c(0,10,20,30,40,50),c(2,2,2,2,2,2),
         lty = 4,col=6)
