# script for BWGS tutorial



#YieldGBLUP <-bwgs.cv (TRAIN47K, YieldBLUE, geno.impute.method="mni", predict.method= "gblup", nFolds=10, nTimes=10 )
#YieldLASSO <-bwgs.cv (TRAIN47K, YieldBLUE, geno.impute.method="mni", predict.method= "LASSO", nFolds=10, nTimes=10 )
#YieldBA <-bwgs.cv (TRAIN47K, YieldBLUE, geno.impute.method="mni", predict.method= "BA", nFolds=10, nTimes=10 )
#YieldRKHS <-bwgs.cv (TRAIN47K, YieldBLUE, geno.impute.method="mni", predict.method= "RKHS", nFolds=10, nTimes=10 )
#YieldEGBLUP <-bwgs.cv (TRAIN47K, YieldBLUE, geno.impute.method="mni", predict.method= "EGBLUP", nFolds=10, nTimes=10 )

#compareM=cbind(YieldGBLUP$cv, YieldLASSO$cv, YieldBA$cv, YieldRKHS$cv, YieldEGBLUP$cv)
#colnames(compareM) = c("GBLUP","LASSO","BayesA","RKHS","EGBLUP")
#boxplot(compareM,xlab="Prediction method",ylab="predictive ability",main="Predictive ability of 5 methods. Yield with 47K markers")


#YieldGBLUP100 <-bwgs.cv (TRAIN47K, YieldBLUE,pop.reduct.method="RANDOM", sample.pop.size=100, geno.impute.method="mni", predict.method="gblup", nFolds=10, nTimes=10 ) 
#YieldGBLUP300 <-bwgs.cv (TRAIN47K, YieldBLUE,pop.reduct.method="RANDOM", sample.pop.size=300, geno.impute.method="mni", predict.method="gblup", nFolds=10, nTimes=10 ) 
#YieldGBLUP500 <-bwgs.cv (TRAIN47K, YieldBLUE, pop.reduct.method="RANDOM",sample.pop.size=500, geno.impute.method="mni", predict.method="gblup", nFolds=10, nTimes=10 )
 
#boxplot(cbind(YieldGBLUP100$cv, YieldGBLUP300$cv, YieldGBLUP500$cv, YieldGBLUP$cv))

#CompareSize=cbind(YieldGBLUP100$cv, YieldGBLUP300$cv, YieldGBLUP500$cv, YieldGBLUP$cv)
#colnames(CompareSize)=c("N=100","N=300","N=500","N=700")
#boxplot(CompareSize,xlab="Training POP size",ylab="Predictive avility",main="Effect of TRAINING POPULATION SIZE")

#testPREDICT_GBLUP=bwgs.predict(geno_train=TRAIN47K,pheno_train=YieldBLUE,geno_target=TARGET47K,MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",r2="NULL",pval="NULL",
#MAP="NULL",geno.impute.method="MNI",predict.method="GBLUP") 

#testPREDICT_EGBLUP=bwgs.predict(geno_train=TRAIN47K,pheno_train=YieldBLUE,geno_target=TARGET47K,MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",r2="NULL",pval="NULL",
#MAP="NULL",geno.impute.method="MNI",predict.method="EGBLUP") 

#testPREDICT_LASSO=bwgs.predict(geno_train=TRAIN47K,pheno_train=YieldBLUE,geno_target=TARGET47K,MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",r2="NULL",pval="NULL",
#MAP="NULL",geno.impute.method="MNI",predict.method="LASSO") 

#testPREDICT_RKHS=bwgs.predict(geno_train=TRAIN47K,pheno_train=YieldBLUE,geno_target=TARGET47K,MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",r2="NULL",pval="NULL",
#MAP="NULL",geno.impute.method="MNI",predict.method="RKHS") 

#testPREDICT_BayesA=bwgs.predict(geno_train=TRAIN47K,pheno_train=YieldBLUE,geno_target=TARGET47K,MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",r2="NULL",pval="NULL",
#MAP="NULL",geno.impute.method="MNI",predict.method="BA") 

#ComparePRED=cbind(testPREDICT_GBLUP[,1] ,testPREDICT_BayesA[,1] ,testPREDICT_LASSO[,1], testPREDICT_RKHS[,1], testPREDICT_EGBLUP[,1])
#colnames(ComparePRED=)c("GBLEP","BauesA","LASSO","RKHS","EGBLUP")
#pairs(ComparePRED,lower.panel = panel.smooth,upper.panel = panel.cor,diag.panel=panel.hist)




TRAIN47K_NO_NA=MNI(TRAIN47K)

datasim03 <- qtlSIM (TRAIN47K_NO_NA, NQTL=100,h2=0.3)
datasim05 <- qtlSIM (TRAIN47K_NO_NA,NQTL=100,h2=0.5)
datasim08 <- qtlSIM(TRAIN47K_NO_NA,NQTL=100,h2=0.8)

 cbind(rownames(datasim03$newSNP),names(datasim03$pheno))





SIM03 <- bwgs.cv (datasim03$newSNP, datasim03$pheno, geno.impute.method="MNI", predict.method ="gblup", nTimes=20, nFolds=5) #
SIM05 <- bwgs.cv (datasim05$newSNP, datasim05$pheno, geno.impute.method="MNI", predict.method="gblup", nTimes=20, nFolds=5) #
SIM08 <- bwgs.cv (datasim08$newSNP, datasim08$pheno, geno.impute.method="MNI", predict.method="gblup", nTimes20, nFolds=5) #
CompareH2=cbind (SIM03,SIM05,SIM08)

colnames(CompareH2)=c("h²=0.3","h²=0.5","h²=0.8")
boxplot(CompareH2,xlab="Simulated Trait heritability",ylab="Predictive avility",main="Effect of TRAIT heritability")


 





