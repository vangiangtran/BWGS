#gstools_pipeline - derived from bwgs
#@ 2017 Sophie Bouchet, Louis Gautier Tran, CHARMET Gilles
#Version 1.0.0 - Release date: 31/02/2017
#bwgs_gilles_v55.R
#sophie.bouchet@inra.fr
#gilles.charmet@inra.fr
#/////////////////////////////////////////////////////////////////////



percentNA=function(x)
{
  perNA=length(x[is.na(x)])/length(x)
  perNA
}

#/////////////////////////////////////////////////////////////////////
#START bwgs.predict TO BE USED ONCE the best model is chosen using bwgs.cv
#/////////////////////////////////////////////////////////////////////

myMean=function(x)
{
  
  DENO=x[!is.na(x)]
  NUME=DENO[DENO>0]
  M=length(NUME)/length(DENO)
  
  M
}





#' Title
#'
#' @param geno_train 
#' @param pheno_train 
#' @param geno_target 
#' @param FIXED_train 
#' @param FIXED_target 
#' @param MAXNA 
#' @param MAF 
#' @param geno.reduct.method 
#' @param reduct.size 
#' @param r2 
#' @param pval 
#' @param MAP 
#' @param geno.impute.method 
#' @param predict.method 
#'
#' @return
#' @export
bwgs.predict <- function(geno_train,pheno_train,geno_target,FIXED_train="NULL",FIXED_target="NULL",MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",r2="NULL",pval="NULL",
                         MAP="NULL",geno.impute.method="NULL",predict.method="GBLUP") 
{
  #(c)2015 louis.gautier.tran@gmail.com & gilles.charmet@clermont.inra.fr
  message("2017 BWGS  - Version 1.6.0 Release date: 31/02/2017")
  #message("2015 Gilles Charmet & Louis Gautier Tran)
  
  # upload the necessary libraries
  
  library(rrBLUP)
  library(BGLR)
  library(glmnet)
  library(e1071) 
  library(randomForest)
  
  
  #message("")
  start.time <- Sys.time()
  message("Start time:")
  print(start.time)
  message("")
  geno.reduct.method <- toupper(geno.reduct.method)
  #reduct.size <- as.numeric(reduct.size)
  geno.impute.method <- toupper(geno.impute.method)
  predict.method <- toupper(predict.method)
  #r2 is for LD shrink method
  #r2 and n.parts are for GMSLD shrink method
  
  #///////////////////////////////////////////////////////////////
  #STEP 0: select common lines in geno and pheno matrices
  #  select common markers in train and target matrices
  # FILTER according to MAF and MAXNA
  #//////////////////////////////////////////////////////////////
  
  if(MAP!="NULL")
  {
    MAPPED_markers=intersect(rownames(MAP),colnames(geno_train))
    
    MAP=MAP[MAPPED_markers,]
    geno_train=geno_train[,MAPPED_markers]
  }
  
  if(FIXED_train!="NULL")
  {
    genoTRFIX=intersect(rownames(FIXED_train),rownames(geno_train))
    
    FIXED_train=FIXED_train[genoTRFIX,]
    geno_train=geno_train[genoTRFIX,]
  }
  
  if(FIXED_target!="NULL")
  {
    genoTAFIX=intersect(rownames(FIXED_target),rownames(geno_target))
    
    FIXED_target=FIXED_target[genoTAFIX,]
    geno_target=geno_target[genoTAFIX,]
  }
  
  
  pheno_train=pheno_train[!is.na(pheno_train)]
  listGENO=rownames(geno_train)
  listPHENO=names(pheno_train)
  LIST=intersect(listGENO,listPHENO)
  
  if (length(LIST)==0)
  {stop("NO COMMON LINE BETWENN geno AND pheno")}
  else
  {
    geno_train=geno_train[LIST,]
    pheno_train=pheno_train[LIST]
    new.pheno.size=length(LIST)
    message("Number of common lines between geno and pheno")
    print(new.pheno.size)
    
    if(FIXED_train!="NULL")
    {
      genoTRFIX=intersect(rownames(FIXED_train),rownames(geno_train))
      FIXED_train=FIXED_train[genoTRFIX,]
      geno_train=geno_train[genoTRFIX,]
    }
    
    
  }
  
  
  marker_train=colnames(geno_train)
  marker_target=colnames(geno_target)
  marker_common=intersect(marker_train,marker_target)
  geno_train=geno_train[,marker_common]
  geno_target=geno_target[,marker_common]
  
  geno <- rbind(geno_train,geno_target)
  geno_train_nrows <- nrow(geno_train)
  geno_target_nrows=nrow(geno_target)
  geno_nrows <- nrow(geno)
  
  if(FIXED_train!="NULL")
  {
    FIXED <- rbind(FIXED_train,FIXED_target)
    FIXED <-MNI(FIXED)
    FIXED <-round(FIXED)
  }
  
  # FILTERING GENOTYPING MATRIX
  # for % NA per marker
  
  markerNA=apply(geno,2,percentNA)
  geno=geno[,markerNA<MAXNA]
  
  freqSNP=apply(geno,2,myMean)
  geno=geno[,freqSNP>MAF & freqSNP<(1-MAF)]
  
  new.geno.size=dim(geno)
  message("Number of markers after filtering")
  print(new.geno.size)
  
  if(MAP!="NULL")
  {
    MAP=MAP[colnames(geno),]
  }
  
  
  
  #//////////////////////////////////////////////////////////////
  #Step 1: Imputation "MNI" or "EMI"
  #//////////////////////////////////////////////////////////////
  
  if((geno.impute.method=="NULL")|(geno.impute.method=="MNI")|(geno.impute.method=="EMI"))
  {
    
    if(geno.impute.method=="MNI"){ 
      
      
      time.mni = system.time(geno_impute <- MNI(geno))
      time.mni.impute = as.numeric(round(time.mni[3]/60,digits=2))
      message("Imputed by MNI.")
      message("Time of imputation by MNI (mins):")
      print(time.mni.impute)
      
      message("A part 5x20 of Imputed genotypic matrix:")
      print(geno_impute[1:5,1:20],quote=FALSE)
      message("")
    } 
    
    
    if(geno.impute.method=="EMI")
    { 
      
      # emi.time.start = system.time()
      # geno_impute <- EMI(geno_shrink)
      # emi.time.stop = system.time()
      time.emi = system.time(geno_impute <- EMI(geno))
      time.emi.impute = as.numeric(round(time.emi[3]/60,digits=2))
      message("Imputed by EMI.")
      message("Time of imputation by EMI (mins):")
      print(time.emi.impute)
      
      message("A part 5x20 of Imputed genotypic matrix:")
      print(geno_impute[1:5,1:20],quote=FALSE)
      message("")
      
    }
    
    message("Imputed by MNI, EMI...finished.")
    geno=geno_impute
    geno_train_impute=geno[1:geno_train_nrows,]
    geno_valid_impute=geno[(geno_train_nrows+1):geno_nrows,]
    
    if(FIXED_train!="NULL")
    {
      FIXED_train=FIXED[1:geno_train_nrows,]
      FIXED_target=FIXED[(geno_train_nrows+1):geno_nrows,]
    } 
  }
  else {stop("Please choose an impute method: NULL, MNI, EMI ")
  } #If not chosen an impute method
  
  
  
  #geno.impute.method can be = "MNI", "EMI"
  #predict.method can be = "LASSO","SVM","BA","BB","BC","RR"(Ridge Regression), "rrBLUP", "RKHS", "MKHS", "EN" (Elastic-Net) 
  #geno.reduct.method can be = "RMR","ANO","LD", "NULL"
  #
  
  #///////////////////////////////////////////////////////////////
  #STEP 2: reduction OF THE DATA (reduire les donnees)
  #//////////////////////////////////////////////////////////////
  
  
  if((geno.reduct.method=="NULL")|(geno.reduct.method=="RMR")|(geno.reduct.method=="ANO")|(geno.reduct.method=="LD")|(geno.reduct.method=="ANO+LD")) 
  { 
    
    if(geno.reduct.method=="NULL"){ 
      geno_shrink <- geno 
      message("No reduction for genomic data.")
      new.geno.size <- dim(geno_shrink)
    }
    
    
    if(geno.reduct.method=="RMR"){   
      
      #if(is.null(reduct.size)) {stop("Please choose the size of columns for new genotypic data.")}
      
      if(!is.null(reduct.size)){
        #shrink_size <- as.numeric(reduct.size)
        time.rmr = system.time(geno_shrink <- RMR(geno,reduct.size))
        time.rmr.reduct = as.numeric(round(time.rmr[3]/60,digits=2))
        new.geno.size <- dim(geno_shrink)
        #geno_shrink
        message("Reduced by RMR. New genotypic data dimension:")
        print(new.geno.size)
        message("")
        message("Time of reduction by RMR (mins):")
        print(time.rmr.reduct)
        
      }  else {stop("Please choose the size of columns for new genotypic data.")}
      
    }
    
    
    if(geno.reduct.method=="ANO"){   
      
      #if(is.null(reduct.size)) {stop("Please choose the size of columns for new genotypic data.")}
      
      if(!is.null(pval)){
        
        time.ano = system.time(genoTrain_shrink <- ANO(pheno_train,geno_train_impute,pval))
        time.ano.reduct = as.numeric(round(time.ano[3]/60,digits=2))
        genoTrain_shrink=genoTrain_shrink[,!is.na(colnames(genoTrain_shrink))]
        genoValid_shrink=geno_valid_impute[,colnames(genoTrain_shrink)]
        geno_shrink=rbind(genoTrain_shrink,genoValid_shrink)
        
        new.geno.size <- dim(geno_shrink)
        
        message("Reduced by ANO. New genotypic data dimension:")
        print(new.geno.size)
        message("")
        message("Time of reduction by ANO (mins):")
        print(time.ano.reduct)
        
      }  else {stop("Please choose the p value for new genotypic data.")}
      
    }
    
    if(geno.reduct.method=="ANO+LD")
    {    
      if(!is.null(pval))
      {
        
        time.ano = system.time(genoTrain_shrinkANO <- ANO(pheno_train,geno_train_impute,pval))
        time.ano.reduct = as.numeric(round(time.ano[3]/60,digits=2))
        genoTrain_shrinkANO=genoTrain_shrinkANO[,!is.na(colnames(geno_shrinkANO))]
        genoValid_shrinkANO=geno_valid_impute[,colnames(genoTrain_shrink)]
        geno_shrinkANO=rbind(genoTrain_shrink,genoValid_shrink)
        
        
      }  else {stop("Please choose the p value for new genotypic data.")}
      
      if(MAP=="NULL")
      {
        stop("Please choose the r2 and/or MAP for LD reduction.")
      } # END OF LD reduction
      
      MAP2=MAP[colnames(geno_shrinkANO),]
      time.chrld = system.time(geno_shrink <- CHROMLD(geno_shrinkANO,R2seuil,MAP2))
      time.chrld.reduct = as.numeric(round(time.chrld[3]/60,digits=2))
      new.geno.size <- dim(geno_shrink)
      #geno_shrink
      time.anold.reduct=time.ano.reduct+time.chrld.reduct
      new.geno.size <- dim(geno_shrink)
      #geno_shrink
      message("Reduced by ANO + CHROMLD. New genotypic data dimension:")
      print(new.geno.size)
      message("")
      message("Time of reduction by ANO+ LD (mins):")
      print(time.anold.reduct)
      
    }
    
    
    if(geno.reduct.method=="LD")
    {
      if(MAP=="NULL")
      {
        stop("Please choose MAP for LD reduction.")
      } # END OF LD reduction
      
      
      MAP2=MAP[colnames(geno_shrinkANO),]
      time.chrld = system.time(geno_shrink <- CHROMLD(geno_impute,R2seuil,MAP))
      time.chrld.reduct = as.numeric(round(time.chrld[3]/60,digits=2))
      new.geno.size <- dim(geno_shrink)
      #geno_shrink
      message("Reduced by CHROMLD. New genotypic data dimension:")
      print(new.geno.size)
      message("")
      message("Time of reduction by CHROMLD (mins):")
      print(time.chrld.reduct)
    }
    
    
    
    # END OF LD reduction
    
    geno=geno_shrink
    geno_train_impute=geno[1:geno_train_nrows,]
    geno_valid_impute=geno[(geno_train_nrows+1):geno_nrows,]
    
    if(FIXED_train!="NULL")
    {
      FIXED_train=FIXED[1:geno_train_nrows,]
      FIXED_target=FIXED[(geno_train_nrows+1):geno_nrows,]
    }
    
  } # End if not chosen a reduct method
  
  else {stop("BWGS Warning! Please choose a valid reduct method. Ex.: RMR, LD, LD+ANO or NULL.")} #If not choose a reduct method
  
  #End of dimension reduction for the genomic matrix 
  
  
  
  
  #/////////////////////////////////////////////////////
  #STEP 3: BREEDING VALUE PREDICTION PROCESS
  #/////////////////////////////////////////////////////
  
  if((predict.method=="EN")|(predict.method=="SVM")|(predict.method=="EGBLUP")|(predict.method=="RR")|(predict.method=="BA")|
     (predict.method=="BB")|(predict.method=="LASSO")|(predict.method=="BL")|(predict.method=="BC")|(predict.method=="GBLUP")|
     (predict.method=="BRNN")|(predict.method=="MKRKHS")|(predict.method=="RF")|(predict.method=="BRR")|(predict.method=="RKHS"))
  {
    
    if (predict.method=="GBLUP") { 
      message("Predict by GBLUP...")
      #GBLUP <-predict_GBLUP(pheno_train,geno_train_impute,geno_valid_impute)
      #GBLUP <-predict_GBLUP(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
      #genoPred=geno_valid_impute, FixedPred=FIXED_target)
      
      GBLUP <-predict_GBLUP(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                            genoPred=geno_valid_impute, FixedPred=FIXED_target)
      
      
      Results <- round(GBLUP,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="EGBLUP") { 
      message("Predict by EGBLUP...")
      EGBLUP <-predict_EGBLUP(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                              genoPred=geno_valid_impute, FixedPred=FIXED_target)
      Results <- round(EGBLUP,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    
    if (predict.method=="BA") { 
      message("Predict by BayesA...")
      BA <-predict_BA(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                      genoPred=geno_valid_impute, FixedPred=FIXED_target)
      Results <- round(BA,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    
    if (predict.method=="BB") { 
      message("Predict by BayesB...")
      BB <-predict_BB(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                      genoPred=geno_valid_impute, FixedPred=FIXED_target)
      Results <- round(BB,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="BC") { 
      message("Predict by BayesC...")
      BC <-predict_BC(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                      genoPred=geno_valid_impute, FixedPred=FIXED_target)
      Results <- round(BC,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="BL") { 
      message("Predict by Bayesian Lasso...")
      BL <-predict_BL(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                      genoPred=geno_valid_impute, FixedPred=FIXED_target)
      Results <- round(BL,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="RF") { 
      message("Predicting by RF...")
      RF <-predict_RF(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                      genoPred=geno_valid_impute, FixedPred=FIXED_target)
      
      message("Predict by RF...ended.")
      Results <- round(RF,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Time stop at:")
      print(end.time)
      print(time.taken,title=FALSE)
      
    } 
    
    if (predict.method=="RKHS") { 
      message("Predicting by RKHS...")
      RKHS <-predict_RKHS(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                          genoPred=geno_valid_impute, FixedPred=FIXED_target)
      Results <- round(RKHS,digits=5)
      #Results = t(Results)
      message("Predict by RKHS...ended.")
      #Results <- round(RKHS,digits=2)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
      
    } 
    
    if (predict.method=="MKRKHS") { 
      message("Predicting by Multi-kernel RKHS...")
      MKRKHS <-
        
        
        predict_MKRKHS(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                       genoPred=geno_valid_impute, FixedPred=FIXED_target)
      Results <- round(MKRKHS,digits=5)
      #Results = t(Results)
      message("Predict by Multi-kernel RKHS...ended.")
      #Results <- round(MKRKHS,digits=2)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Stop time:")
      print(end.time)
      print(time.taken,title=FALSE)
      
    } 
    
    
    if (predict.method=="RR") { 
      message("Predicting by RR...")
      RR <-predict_RR(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                      genoPred=geno_valid_impute, FixedPred=FIXED_target)
      # RR = t(RR) 
      #rownames(RR)="";
      message("Predict by RR...ended.")
      Results <- round(RR,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Time stop at:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="BRR") { 
      message("Predicting by BRR...")
      BRR <-predict_BRR(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                        genoPred=geno_valid_impute, FixedPred=FIXED_target)
      #BRR = t(BRR) 
      #rownames(RR)="";
      message("Predict by BRR...ended.")
      Results <- round(BRR,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Time stop at:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    
    if (predict.method=="BRNN") { 
      message("Predicting by BRNN...")
      BRNN <-predict_BRNN(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                          genoPred=geno_valid_impute, FixedPred=FIXED_target)
      
      message("Predict by BRNN...ended.")
      Results <- round(BRNN,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Time stop at:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="LASSO") { 
      message("Predicting by LASSO...")
      LASSO <-predict_Lasso(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                            genoPred=geno_valid_impute, FixedPred=FIXED_target)
      #LASSO = t(LASSO) 
      #rownames(LASSO)="";
      #rownames(LASSO)=rownames(geno_valid_impute)
      message("Predict by LASSO...ended.")
      Results <- round(LASSO,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Time stop at:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="EN") { 
      message("Predicting by Elastic-Net...")
      ElasticNet <-predict_ElasticNet(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                                      genoPred=geno_valid_impute, FixedPred=FIXED_target)
      #ElasticNet = t(ElasticNet)
      #rownames(ElasticNet)=""
      message("Predict by Elastic-Net...ended.")
      Results <- round(ElasticNet,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Time stop at:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    if (predict.method=="SVM") { 
      message("Predicting by SVM...")
      SVM <-predict_SVM(phenoTrain=pheno_train, genoTrain=geno_train_impute, FixedTrain=FIXED_train, 
                        genoPred=geno_valid_impute, FixedPred=FIXED_target)
      #SVM = t(SVM)
      #rownames(SVM)=""
      message("Predict by SVM...ended.")
      Results <- round(SVM,digits=5)
      #rownames(Results) <- "";
      message("Phenotypic estimation GEBV:")
      message("")
      print(Results)
      message("")
      
      end.time <- Sys.time()
      time.taken <- end.time-start.time
      message("Time stop at:")
      print(end.time)
      print(time.taken,title=FALSE)
    } 
    
    
    #  } #else stop("Please choose an impute method: MNI, EMI")}
    
    # else {stop("Please choose a predict method: EN, SVM, RR, LASSO, GBLUP, RKHS, BA, BB, BC, BL, RF")}
    
    
  } else {
    stop("Please choose a predict method: EN, SVM, BRNN, BRR, RR, LASSO, GBLUP,EGBLUP, RKHS, MKRKHS,BA, BB, BC, BL, RF")
  } # End of all predict methods
  
  return(Results)
  
} #END OF BWGS.PREDICT()




##REDUCTION MODELS and other miscellaneous functions

#////////////////////////////////////////////////////////////
#Imputation by the A matrix Jeffrey Endelman JF.Shrink()
#////////////////////////////////////////////////////////////

AM <- function(geno,arg.A.mat=list(min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,n.core=1,shrink=TRUE,return.imputed=TRUE))
{
  
  
  #Shrink_result = A.mat(X=geno, impute.method = "mean", tol = 0.02, return.imputed = FALSE)$A
  #library(rrBLUP)
  Shrink_result = A.mat(X=geno, impute.method = "mean", tol = 0.02,n.core=1,shrink=TRUE,return.imputed = TRUE)$A
  # Shrink_result = A.mat,args=c(X=geno,arg.A.mat))$A
  return(Shrink_result)  
}


#///////////////////////////////////////
#  RMR: Randaom sampling of markers
#///////////////////////////////////////

RMR <- function(geno,N) 
  
{
  ncolG <- ncol(geno)
  V <- sample(1:ncolG,N,replace=F)
  random_geno <- geno[,V]
  return(random_geno)
}

#//////////////////////////////////////////////////////////////
#   ANO: selection of N markers with the lowest Pvalue in GWAS
#//////////////////////////////////////////////////////////////

ANO <- function(P,GG,pval) 
{
  
  #Genotypic Reduction by Association (ANOVA)
  #(c)15/06/2015 G. CHARMET & V.G. TRAN
  #GG = geno
  #P = pheno
  #genotypes and phenotypes must not have NAs
  
  listF = ANOV(P,GG)
  names(listF)=colnames(GG)
  ChoixM = GG[,listF < pval&!is.na(listF)]
  dim(ChoixM)
  ChoixM
}

# function to run ANOVA
ANOV <- function(P,GG)
{
  TESTF=rep(0,ncol(GG))
  for (i in 1:ncol(GG))
  {
    GGi=GG[,i]
    names(GGi)=rownames(GG)
    GGii=GGi[!is.na(GGi)]
    Pi=P[!is.na(GGi)]
    AOV=anova(lm(Pi~GGii))
    TEST=AOV["Pr(>F)"][1]
    TESTF[i]=c(TEST)[[1]][1]
  }
  return(TESTF)
}

#////////////////////////////////////////////////////////////////////////////////////
# LD  function to select markers based on LD: eliminate pairs with highest LD rd
# NB: does NOT work: eliminate TOO MANY markers
# MOREOVER LDCORSV is TOO slow
#////////////////////////////////////////////////////////////////////////////////////



#////////////////////////////////////////////////////////////////////////////////////////
#
#  NEWLD is a faster function to remove markers in LD > lambda
#
#/////////////////////////////////////////////////////////////////////////////////////////

NEWLD <- function(geno,R2seuil)
{
  genoNA=MNI(geno)
  CORLD=cor(genoNA,use="pairwise.complete.obs")
  CORLD=abs(CORLD)
  diag(CORLD)=0
  CORLD[is.na(CORLD)]<-0
  #	CORLD[upper.tri(CORLD)]<-0
  range(CORLD)
  geno_LD=CORLD
  R2maxRow=apply(geno_LD,1,max)
  R2maxCol=apply(geno_LD,2,max)
  
  geno_LD=geno_LD[R2maxRow<R2seuil,R2maxCol<R2seuil]
  #geno_final=NEWLD_clean(geno,CORLD,R2seuil=R2seuil)
  geno_LD=geno[,colnames(geno_LD)]
  
  return(geno_LD)
}

#/////////////////////////////////////////////////////////////////////////////////////////////
#NEWLD clean function: another function for stepwise elimination of markers with the highest LD
#  can be VERY llong for large matrices
#/////////////////////////////////////////////////////////////////////////////////////////////

# @param geno_data is the genotyping matrix with markers in columns
#

NEWLD_clean <- function(geno_data,CORLD,R2seuil) 
  
{
  #  R2seuil=lambda
  # function to select a subset of marker from columns of geno_data
  # so that maximum r2 between pairs of markers is < R2seuil
  # LD is output matrix from LD.Measures of package LDcorSV with supinfo=TRUE
  # typeR2 is the type of R2: "N" for "normal", "V" for relationship corrected 
  # "S" for structure corrected and "SV" for both corrections S+V
  
  marker2remove=character()# to cumulate markers to be removed from the geno_data
  
  geno_LD=CORLD
  # R2max=max(geno_LD)
  
  R2maxRow=apply(geno_LD,1,max)
  R2maxCol=apply(geno_LD,2,max)
  
  geno_LD=geno_LD[R2maxRow<R2seuil,R2maxCol<R2seuil]
  
  
  
  #geno_LD
  
  Newgeno_data<-geno_data[,colnames(geno_LD)]
  
  #output
  return(Newgeno_data)
}






#/////////////////////////////////////////////////////////////////////////////
# CHROMLD: function to split a genotypic matrix into chromosomes
# then apply LD within each chrmosome (save time)
#/////////////////////////////////////////////////////////////////////////////

CHROMLD <- function(geno,R2seuil,MAP) 
{ 
  geno_test <- geno
  chrom=MAP[,"chrom"]
  
  Nchrom <- unique(chrom)
  
  M = list()
  E = list()
  
  for (i in Nchrom)
  {
    M[[i]] = geno_test[,chrom==i]
    E[[i]] = NEWLD(M[[i]],R2seuil)
    
  }
  
  E_tout <- do.call(cbind,E) # See function: do.call(rbind,"list of matrices")
  
  #E_tout <- LD(E_tout,lambda) # LD them mot lan nua cua matran ghep de loai bo cac cap marker co LD trong N matrices
  
  return(E_tout)
  
}
#

#////////////////////////////////////////////////////////////////////////////////////
#  RMGG  to create random NA value in geno, useful for testing Imputation methdos
#////////////////////////////////////////////////////////////////////////////////////

RMGG <-function(G, N) 
{
  #RMMG - Random Missing Matrix Generator
  #(c)2014 V.G. TRAN & G. CHARMET
  #G: genotypic matrix
  #N: % missing
  nRow = dim(G)[1]
  nCol = dim(G)[2]
  N_missing = N*nRow*nCol/100
  Gm<-G
  Missing<-sample(nRow*nCol, N_missing, replace=FALSE)
  ROWS<-matrix(data = Gm, ncol=1)
  ROWS[Missing]=NA
  Gm<-matrix(data = ROWS, nrow=nRow, ncol=nCol, byrow=FALSE)
  rownames(Gm) <- rownames(G)
  colnames(Gm) <- colnames(G)
  return(Gm)
}


#////////////////////////////////////////////////////////
#  RPS  function for randomly sampling individuals (lines)
# ONLY useful for teachning purpose, NOT in real life
#////////////////////////////////////////////////////////////


RPS <- function(geno,N) 
{ 
  #Random Pop Size
  #(c)03/12/2014 vangiang.tran@gmail.com & gilles.charmet@clermont.inra.fr
  nrowG <- nrow(geno)
  V <- sample(1:nrowG,N,replace=F)
  random_geno <- geno[V,]
  return(random_geno)
}


#//////////////////////////////////////////////////////////////////////////
#///IMPUTATION MODELS
#/////////////////////////////////////////////////////////////////////////

#/////////////////////////////////////////////////////////////////////////////////////
# EMI: imputation by Expactation-Maximization algorithm
# uses A.mat function from rrBLUP package
#/////////////////////////////////////////////////////////////////////////////////////

#EMI <- function(geno,arg.A.mat=list(impute.method = "EM", tol = 0.02, return.imputed = TRUE)) 
EMI <- function(geno) 
{
  library(rrBLUP)
  
  
  #EMI_result= do.call(rrBLUP::A.mat,args=c(X=geno, arg.A.mat))$imputed
  EMI_result= A.mat(X=geno, impute.method = "EM", tol = 0.02, return.imputed = TRUE)$imputed
  return(EMI_result)
  
}



#////////////////////////////////////////////////////////////////
#CROSS VALIDATION PROCESSING
#///////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////
# runCrossVal  carries out cross validation, return only predicted value, to be run with tha "ALL" option
#//////////////////////////////////////////////////////////////////////////////////////

runCrossVal <- function(pheno, geno, predictor, nFolds, nTimes) 
{
  
  #start.time.cv <- Sys.time()
  
  listGENO=rownames(geno)
  listPHENO=names(pheno)
  LIST=intersect(listGENO,listPHENO)
  
  if (length(LIST)==0)
  {stop("NO COMMON LINE BETWENN geno AND pheno")}
  else
  {
    geno=geno[LIST,]
    pheno=pheno[LIST]
    new.pheno.size=length(LIST)
    message("Number of common lines between geno and pheno")
    print(new.pheno.size)
  }
  
  notIsNA <- !is.na(pheno) # 
  pheno <- pheno[notIsNA] ## 
  geno <- geno[notIsNA,] ## 
  nObs <- length(pheno)
  #saveAcc <- rep(NA, nTimes)
  saveAcc <- rep(NA, nTimes*nFolds)
  savePred <- rep(0, nObs)
  
  # n_row = nTimes*nFolds
  # n_col = 3
  # cvtable <- matrix(rep(0,n_row*n_col),n_row,n_col)
  ###  title_table <- c("Time","Fold","CV correlation");
  ACCU <- rep(NA,nFolds) #initial de ACC pour plusieurs folds 04/05/2015 VG TRAN
  
  for (time in 1:nTimes) {
    
    #   timebar <- time*100/nTimes # For Windows Time Bar
    #   waitingbar(runif(timebar)) # Windows Bar
    
    cvPred <- rep(NA, nObs)
    folds <- sample(rep(1:nFolds, length.out=nObs))
    
    for (fold in 1:nFolds) {
      phenoTrain <- pheno[folds != fold]
      genoTrain <- geno[folds != fold,]
      genoPred <- geno[folds == fold,]
      pred <- predictor(phenoTrain, genoTrain, genoPred)
      #message(pred) # them vao 21/11/14
      cvPred[folds == fold] <- pred
      valid <- pheno[folds == fold] # pheno valid them vao 21/11/14
      # Correlation entre valid/pred 
      #accu <- round(cor(pred,valid),digits=2) # them vao 21/11/2014
      corr <- cor(pred,valid)
      ACCU[fold] <- round(corr[1],digits=2)# 04/05/2015 VG TRAN
      
      message("")  # them vao 21/11/2014
      #                 # cat(c("cv for fold",fold,"is:",accu)) # them vao 21/11/2014
      #message(c("cv correlation for fold ",fold," is: ",accu)) # them vao 21/11/2014
      message(c("cv correlation for fold ",fold," is: ",ACCU[fold])) # 04/05/2015 VG TRAN
      
      #           Save to a matrix time, fold, and cv correlation: 
      #  
      ##      cvTableLine <- c(time,fold,accu)
      ##     cvtable <- rbind(title_table,cvTableLine)
      #names(cvtable) = c("Time","Fold","CV correlation")
      #print(cvtable) lier "cvtable" avec "Results"
      #print(cvtable)
      # SD <- round(sd(accu),digits=4)
    }
    timebar <- time*100/nTimes # For Windows Time Bar
    waitingbar(runif(timebar)) # Windows Bar
    #saveAcc[time] <- round(cor(cvPred, pheno)*100,digits=2)  # pourcentage/100%
    #saveAcc[time] <- round(cor(cvPred, pheno),digits=2) # normal
    #saveAcc[time] <- round(mean(ACCU),digits=3) # 04/05/2015 VG TRAN
    saveAcc[(((time-1)*nFolds)+1):(time*nFolds)] <- round(ACCU,digits=2) # 04/05/2015 VG TRAN
    message(c("Mean CV correlation for time ",time," and ",nFolds," folds is: ",saveAcc[time]))
    savePred <- savePred + cvPred 
    # SD <- round(sd(saveAcc[time]),digits=4)
  } 
  SD <- round(sd(saveAcc),digits=4) #Standard Deviation
  predMean <- rep(NA, length(notIsNA))
  predMean[notIsNA] <- savePred/nTimes
  message("")
  
  #end.time.cv <- Sys.time()
  #time.taken.cv <- end.time.cv-start.time.cv
  #return(c(saveAcc,SD,cvtable)) # 
  return(c(saveAcc,SD)) # OK
  #return(c(cvtable))
  # return(c(saveAcc,SD,time.taken.cv)) # OK
  #Nota: saveAcc: correlation between pheno pred and pheno real, SD: Standard Deviation(?cartype)
}



#////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# runCV  carries out cross validation, returns prediction and CD for each line (available only with some methods)
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

runCV<- function(pheno, geno, FIXED, pop.reduct.method,rps, predictor, nFolds, nTimes) 
{
  #(c)2014 V.G. TRAN & D. LY
  class_cv <- list()#start.time.cv <- Sys.time()
  notIsNA <- !is.na(pheno)
  pheno <- pheno[notIsNA]
  geno <- geno[notIsNA,]
  nObs <- length(pheno)
  saveAcc <- rep(NA, nTimes)
  saveMSEP <- rep(NA, nTimes)
  #saveAcc <- rep(NA, nTimes*nFolds)
  savePred <- rep(0, nObs)
  savePredSD <-rep(0,nObs)
  saveCD <-rep (0,nObs)
  # cvtable <- matrix(rep(0,n_row*n_col),n_row,n_col)
  # ACCU <- rep(0,nFolds) #initial de ACC pour plusieurs folds 04/05/2015 VG TRAN
  
  
  class_predict = list()
  class_predSD = list()
  class_CD = list()
  
  # Predict <- rep(0, nObs)
  #CD <- rep(0,nObs) #Coefficient of Determination
  
  # names(Predict) <- names(pheno)
  #names(CD) <- names(pheno)
  
  for (time in 1:nTimes) 
  {
    #   timebar <- time*100/nTimes # For Windows Time Bar#   waitingbar(runif(timebar)) # Windows Bar#   
    #Random Pop Size (added on 08/06/2015 V.G. TRAN)
    ACCU <- rep(0,nFolds) #initial de ACC pour plusieurs folds 04/05/2015 VG TRAN
    
    
    if (pop.reduct.method=="NULL") 
    {
      # nObs <- length(pheno)
      message("Random Pop Size  is not applied.")
      
      Predict <- rep(0, nObs)
      PredSD <-rep(0, nObs)
      CD <-rep(0, nObs)
      
      
      
      
      folds <- sample(rep(1:nFolds, length.out=nObs))
      
      for (fold in 1:nFolds)
      {
        
        phenoTrain <- pheno[folds != fold]
        phenoReal <- pheno[folds == fold] 
        genoTrain <- geno[folds != fold,]
        genoPred <- geno[folds == fold,]
        if (FIXED!="NULL")
        {
          FixedPred <- FIXED[folds == fold,]
          FixedTrain <- FIXED[folds!=fold,]
          pred <- predictor(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
        }
        
        # pred <- predictor(phenoTrain,phenoReal, genoTrain, genoPred)
        else {pred <- predictor(phenoTrain, genoTrain, FixedTrain="NULL", genoPred, FixedPred="NULL")}
        
        
        Predict[folds == fold] <- pred[,1]
        PredSD[folds==fold] <- pred[,2]
        CD[folds==fold] <- pred[,3]
        
        valid <- pheno[folds == fold] # pheno valid them vao 21/11/14 valid = phenoReal
        # Correlation entre valid/pred 
        
        
        corr <- cor(pred[,1],valid)
        ACCU[fold] <- round(corr[1],digits=3)# 04/05/2015 VG TRAN
        
        message("")
        # message(c("cv correlation for fold ",fold," is: ",accu));
        message(c("cv correlation for fold ",fold," is: ",ACCU[fold])) # 04/05/2015 VG TRAN
        
        # Predict[rownames(pred)] <- pred[,1]
        # PredSD [rownames(pred)] <- pred[,2]
        # CD[rownames(pred)] <-pred[,3]
        
      }  # end of fold loop
      
      message("")
      class_predict[[time]] <- Predict
      class_predSD[[time]] <- PredSD
      class_CD[[time]] <-CD
      timebar <- time*100/nTimes
      waitingbar(runif(timebar))
      
      #saveAcc[time] <- round(cor(cvPred, pheno)*100,digits=2);  # pourcentage/100%
      saveAcc[time] <- round(cor(Predict, pheno, use="na.or.complete"),digits=3) # normal
      saveMSEP[time] <- round(sqrt(((Predict-pheno)^2)/length(pheno)),3)
      # saveAcc[(((time-1)*nFolds)+1):(time*nFolds)] <- round(ACCU,digits=2) # 04/05/2015 VG TRAN
      #saveAcc[time] <- round(mean(ACCU,na.rm=T),digits=2) 
      message("")
      message(c("Mean CV correlation for time ",time," and ",nFolds," folds is: ",saveAcc[time]))
      savePred <- savePred + Predict
      savePredSD <-savePredSD + PredSD
      saveCD <- saveCD+CD
      #savePpred=Predict
      #savePredSD=PredSD
      #saveCD=CD
      
      rm(Predict,PredSD,CD)
    } 
    
    
    
    if(pop.reduct.method=="RANDOM")  
    {              
      genopop <- RPS(geno,rps) # Random Pop Size #ATTENTION: car le size de genopop est different que rps!! POURQUOI
      new.geno.pop.size <- dim(genopop)
      phenopop <- pheno[rownames(genopop)]
      nObs <- rps
      
      message("Random Pop Size applied. New genotypic data dimension:")
      print(new.geno.pop.size)
      message("")
      
      nObs <- length(pheno)
      
      
      Predict<-rep(NA,nObs)
      CD<-rep(NA,nObs)
      PredSD <-rep(NA,nObs)
      
      names(PredSD)=names(pheno)
      names(CD)=names(pheno)
      names(Predict)=names(pheno)
      
      folds <- sample(rep(1:nFolds, length.out=nObs))
      
      for (fold in 1:nFolds)
      {
        
        genoTrain <- RPS(geno,rps)
        phenoTrain=pheno[rownames(genoTrain)]
        phenoReal <- pheno[!names(pheno)%in%names(phenoTrain)] 
        genoPred <- geno[names(phenoReal),]
        if (FIXED!="NULL")
        {
          FixedPred <- FIXED[folds == fold,]
          FixedTrain <- FIXED[folds!=fold,]
          pred <- predictor(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
        }
        
        
        
        #pred <- predictor(phenoTrain,phenoReal, genoTrain, genoPred)
        else {pred <- predictor(phenoTrain, genoTrain, FixedTrain="NULL", genoPred, FixedPred="NULL")}
        #else {pred <- predictor(phenoTrain, genoTrain, genoPred)}
        #Return cbind(GPRED,CD)
        #pred <- predictor_rrBLUP(phenoTrain,phenoReal, genoTrain, genoPred); 
        #Return cbind(GPRED,CD)
        
        
        Predict[rownames(pred)] <- pred[,1]
        
        PredSD[rownames(pred)] <-pred[,2]
        CD[rownames(pred)] <- pred[,3]
        
        valid <- pheno[rownames(pred)] # pheno valid them vao 21/11/14 valid = phenoReal
        # Correlation entre valid/pred 
        
        # accu <- round(cor(pred[,1],phenoReal),digits=2);
        # ACCU[fold] <- round(cor(pred,valid),digits=2)# 04/05/2015 VG TRAN
        corr <- cor(pred[,1],valid)
        ACCU[fold] <- round(corr[1],digits=3)# 04/05/2015 VG TRAN
        
        message("")
        # message(c("cv correlation for fold ",fold," is: ",accu));
        message(c("cv correlation for fold ",fold," is: ",ACCU[fold])) # 04/05/2015 VG TRAN
        
        
      }  # fold loop
      message("")
      class_predict[[time]] <- Predict
      class_predSD[[time]] <- PredSD
      class_CD[[time]] <-CD
      timebar <- time*100/nTimes
      waitingbar(runif(timebar))
      
      #saveAcc[time] <- round(cor(cvPred, pheno)*100,digits=2);  # pourcentage/100%
      saveAcc[time] <- round(cor(Predict, pheno,use="na.or.complete"),digits=3) # normal
      saveMSEP[time] <- round(sqrt(((Predict-pheno)^2)/length(pheno)),3)
      #saveAcc[(((time-1)*nFolds)+1):(time*nFolds)] <- round(ACCU,digits=2) # 04/05/2015 VG TRAN
      #saveAcc[time] <- round(mean(ACCU,na.rm=T),digits=2) 
      message("")
      message(c("Mean CV correlation for time ",time," and ",nFolds," folds is: ",saveAcc[time]))
      savePred <- savePred + Predict
      savePredSD <- savePredSD + PredSD
      saveCD <- saveCD+CD	   
      
      rm(Predict,PredSD,CD)	   
      
      
    }   #//////End of Random Pop Size script.
    
    
    if(pop.reduct.method=="OPTI")  
    {              
      genopop <- optiTRAIN(geno,rps,Nopti=3000)$genoOptimiz # optimum Pop Size #ATTENTION: car le size de genopop est different que rps!! POURQUOI
      new.geno.pop.size <- c(rps,ncol(geno))
      phenopop <- pheno[rownames(genopop)]
      nObs <- rps
      
      message("CDmean optimisation applied. New genotypic data dimension:")
      print(new.geno.pop.size)
      message("")
      
      nObs <- length(pheno)
      
      
      #Predict <- rep(0,length(pheno))
      Predict<-rep(NA,nObs)
      CD<-rep(NA,nObs)
      PredSD <-rep(NA,nObs)
      #cvPred <- rep(NA, nObs)
      #cvPredSD <-rep(NA,nObs)
      #cvCD <- rep(NA,nObs)
      #names(Pred)=names(pheno)
      names(PredSD)=names(pheno)
      names(CD)=names(pheno)
      names(Predict)=names(pheno)
      folds <- sample(rep(1:nFolds, length.out=nObs))
      
      for (fold in 1:nFolds)
      {
        
        testOPTI<- optiTRAIN(geno,rps,Nopti=30) 
        genoTrain <- testOPTI$genoOptimiz 
        phenoTrain=pheno[rownames(genoTrain)]
        phenoReal <- pheno[!names(pheno)%in%names(phenoTrain)] 
        genoPred <- geno[names(phenoReal),]
        
        if (FIXED!="NULL")
        {
          FixedPred <- FIXED[rownames(genoTrain),]
          FixedTrain <- FIXED[!rownames(FIXED)%in%names(phenoTrain),]
          pred <- predictor(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
        }
        
        
        
        #pred <- predictor(phenoTrain,phenoReal, genoTrain, genoPred)
        else {pred <- predictor(phenoTrain, genoTrain, FixedTrain="NULL", genoPred, FixedPred="NULL")}
        
        
        #pred <- predictor(phenoTrain, genoTrain, genoPred)
        
        #Return cbind(GPRED,CD)
        #pred <- predictor_rrBLUP(phenoTrain,phenoReal, genoTrain, genoPred); 
        #Return cbind(GPRED,CD)
        
        
        Predict[rownames(pred)] <- pred[,1]
        PredSD[rownames(pred)] <- pred[,2]
        CD[rownames(pred)] <- pred[,3]
        
        valid <- pheno[rownames(pred)] # pheno valid them vao 21/11/14 valid = phenoReal
        # Correlation entre valid/pred 
        
        # accu <- round(cor(pred[,1],phenoReal),digits=2);
        # ACCU[fold] <- round(cor(pred,valid),digits=2)# 04/05/2015 VG TRAN
        #corr <- cor(pred,valid)
        corr <-cor(pred[,1],valid)
        ACCU[fold] <- round(corr[1],digits=3)# 04/05/2015 VG TRAN
        
        message("")
        # message(c("cv correlation for fold ",fold," is: ",accu));
        message(c("cv correlation for fold ",fold," is: ",ACCU[fold])) # 04/05/2015 VG TRAN
        
        
        
      }
      message("")
      class_predict[[time]] <- Predict
      class_predSD[[time]] <-PredSD
      class_CD[[time]] <-CD
      timebar <- time*100/nTimes
      waitingbar(runif(timebar))
      
      
      
      
      #saveAcc[time] <- round(cor(cvPred, pheno)*100,digits=2);  # pourcentage/100%
      saveAcc[time] <- round(cor(cvPred, pheno,use="pairwise.complete.obs"),digits=3) # normal
      
      # saveAcc[(((time-1)*nFolds)+1):(time*nFolds)] <- round(ACCU,digits=2) # 04/05/2015 VG TRAN
      #saveAcc[time] <- round(mean(ACCU,na.rm=T),digits=3)
      #saveAcc[time] <- round(cor(Predict, pheno[names(Predict)],use="na.or.complete"),digits=3)
      saveMSEP[time] <- round(sqrt(sum((Predict-pheno)^2)/length(pheno)),3)
      message("")
      message(c("Mean CV correlation for time ",time," and ",nFolds," folds is: ",saveAcc[time]))
      savePred <- savePred + Predict
      savePredSD <- savePredSD + PredSD
      saveCD <- saveCD+CD	  
      rm(Predict,PredSD,CD)    
      
    }   #//////End of OPTI script.
    
  }# end nTimes
  
  SD <- round(sd(saveAcc),digits=4) # Standard Deviation of CV
  SDMSEP=round(sd(saveMSEP),4)      # Standard Deviation of sqrMSEP
  
  predMean <- rep(NA, length(notIsNA))
  predMean[notIsNA] <- savePred / nTimes
  
  CDMean <- rep(NA, length(notIsNA))
  CDMean[notIsNA] <- saveCD / nTimes
  
  
  
  message("")
  #  message("CV correlations for all times : ")
  #  print(saveAcc);
  #bv_table_mean <- mean(CLASS_TABLE)
  #message("BV table mean : ")     
  
  #print(class_predict)
  
  bv_predict_all = do.call(cbind,class_predict)
  bv_predSD_all = do.call(cbind,class_predSD)
  bv_CD_all = do.call(cbind,class_CD)
  # message("All BV predictions: ")
  # print(bv_predict_all)
  bv_predict_mean = apply(bv_predict_all,1,meanNA)
  bv_predSD_mean=apply(bv_predSD_all,1,meanNA)
  # bv_predict_mean=bv_predict_mean+mean(pheno)
  bv_CD_mean= apply(bv_CD_all,1,meanNA)
  
  #message("BV predict mean: ");
  #print(bv_predict_mean);
  #message("BV Predict all: ");      
  #print(bv_predict_all)
  #message("BV Predict mean: "); 
  #print(bv_predict_mean);     
  
  #mae <- abs(bv_predict_mean-pheno)*100/pheno
  bv_table <- round(cbind(pheno,bv_predict_mean,bv_predSD_mean,bv_CD_mean),digits=3)
  #rownames(bv_table) <- names(pheno)
  #colnames(bv_table) <- c("Real pheno","Predict BV","predict BV_SD","CDmean")
  class_cv[[1]] <- c(saveAcc,SD)
  class_cv[[2]] <- c(saveMSEP,SDMSEP)
  class_cv[[3]] <- bv_table
  
  
  return(class_cv)
  
  #return(class_cv)
  
}
meanNA <- function(x)
{ 
  m=mean(x[!is.na(x)])
  m
}

#//////////////////////////////////////////////////////////////////////////////////
#  Run_All_Croos_Validation
#  performs prediction with a range of methods (not all)
#  and provide a summary table for comparison
#/////////////////////////////////////////////////////////////////////////////////////

Run_All_Cross_Validation <- function(pheno,geno_impute,nFolds,nTimes)
{
  #(c)2013 V.G. TRAN
  
  #cross validation for predict_ElasticNet(with glmnet): 
  
  #Elastic-Net sequence (alpha not equal to 1 and zero.
  message("Predicting by Elastic-Net ...")
  time.en = system.time(ElasticNet <-runCrossVal(pheno,geno_impute,predict_ElasticNet,nFolds,nTimes))
  L = length(ElasticNet)
  message("Predict by Elastic-Net...ended.")
  EN_Results <- c(mean(ElasticNet[1:L-1],na.rm=T),ElasticNet[L],round(time.en[3]/60,digits=3))
  names(EN_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(EN_Results)
  message("")
  
  
  #cross validation for predictor_SVM:
  
  message("Predicting by SVM ...")
  time.en = system.time(SVM <-runCrossVal(pheno,geno_impute,predict_SVM,nFolds,nTimes))
  L = length(ElasticNet)
  message("Predict by SVM...ended.")
  SVM_Results <- c(mean(SVM[1:L-1],na.rm=T),SVM[L],round(time.en[3]/60,digits=3))
  names(SVM_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(SVM_Results)
  message("")
  
  
  #cross validation for predictor_RR(Ridge Regression penalty with glmnet): 
  #Ridge Regression Model
  
  message("Predicting by RR (Ridge Regression)...")
  time.rr = system.time(RR <-runCrossVal(pheno,geno_impute,predict_RR,nFolds,nTimes))
  L = length(RR)
  message("Predict by RR...ended.")
  RR_Results <- c(mean(RR[1:L-1],na.rm=T),RR[L],round(time.rr[3]/60,digits=3))
  names(RR_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(RR_Results)
  message("")
  
  
  #cross validation for predictor_Lasso(LASSO penalty with glmnet):
  
  message("Predicting by LASSO ...")
  time.lasso = system.time(LASSO <-runCrossVal(pheno,geno_impute,predict_Lasso,nFolds,nTimes))
  L = length(LASSO)
  message("Predict by LASSO...ended.")
  LASSO_Results <- c(mean(LASSO[1:L-1],na.rm=T),LASSO[L],round(time.lasso[3]/60,digits=3))
  names(LASSO_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(LASSO_Results)
  message("")
  
  
  #cross validation for predictor_GBLUP:
  
  message("Predict by GBLUP...")
  time.gblup = system.time(GBLUP <-runCrossVal(pheno,geno_impute,predict_GBLUP,nFolds,nTimes))
  L = length(GBLUP)
  message("Predict by GBLUP...ended.")
  GBLUP_Results <- c(mean(GBLUP[1:L-1],na.rm=T),GBLUP[L],round(time.gblup[3]/60,digits=3))
  names(GBLUP_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(GBLUP_Results)
  message("")
  
  #cross validation for predictor_EGBLUP:
  
  message("Predict by EGBLUP...")
  time.egblup = system.time(EGBLUP <-runCrossVal(pheno,geno_impute,predict_EGBLUP,nFolds,nTimes))
  L = length(EGBLUP)
  message("Predict by EGBLUP...ended.")
  EGBLUP_Results <- c(mean(EGBLUP[1:L-1],na.rm=T),EGBLUP[L],round(time.egblup[3]/60,digits=3))
  names(EGBLUP_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(EGBLUP_Results)
  message("")
  
  
  #cross validation for predictor_BA (Bayes A):
  message("Predicting by BA...")
  time.BA = system.time(BA <-runCrossVal(pheno,geno_impute,predict_BA,nFolds,nTimes))
  L = length(BA)
  message("Predict by Bayes A...ended.")
  BA_Results <- c(mean(BA[1:L-1],na.rm=T),BA[L],round(time.BA[3]/60,digits=3))
  names(BA_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(BA_Results)
  message("")
  
  #cross validation for predictor_BB (Bayes B):
  message("Predicting by BB...")
  time.BB = system.time(BB <-runCrossVal(pheno,geno_impute,predict_BB,nFolds,nTimes))
  L = length(BB)
  message("Predict by Bayes B...ended.")
  BB_Results <- c(mean(BB[1:L-1],na.rm=T),BB[L],round(time.BB[3]/60,digits=3))
  names(BB_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(BB_Results)
  message("")
  
  #cross validation for predictor_BC (Bayes C):
  message("Predicting by BC...")
  time.BC = system.time(BC <-runCrossVal(pheno,geno_impute,predict_BC,nFolds,nTimes))
  L = length(BC)
  message("Predict by Bayes C...ended.")
  BC_Results <- c(mean(BC[1:L-1],na.rm=T),BC[L],round(time.BC[3]/60,digits=3))
  names(BC_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(BC_Results)
  message("")
  
  #cross validation for predictor_BL (Bayesian Lasso):
  message("Predicting by BL...")
  time.BL = system.time(BL <-runCrossVal(pheno,geno_impute,predict_BL,nFolds,nTimes))
  L = length(BL)
  message("Predict by Bayes LASSO...ended.")
  BL_Results <- c(mean(BL[1:L-1],na.rm=T),BL[L],round(time.BL[3]/60,digits=3))
  names(BL_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(BL_Results)
  message("")
  
  #cross validation for predictor_BRR (Bayes RR):
  message("Predicting by BRR...")
  time.BRR = system.time(BRR <-runCrossVal(pheno,geno_impute,predict_BRR,nFolds,nTimes))
  L = length(BRR)
  message("Predict by Bayes RR...ended.")
  BRR_Results <- c(mean(BRR[1:L-1],na.rm=T),BRR[L],round(time.BRR[3]/60,digits=3))
  names(BRR_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(BRR_Results)
  message("")
  
  #cross validation for predictor_BRNN (Bayesian regularization Neural Network):
  message("Predicting by BRNN...")
  time.BRNN = system.time(BRNN <-runCrossVal(pheno,geno_impute,predict_BRNN,nFolds,nTimes))
  L = length(BRNN)
  message("Predict by Bayes RR...ended.")
  BRNN_Results <- c(mean(BRNN[1:L-1],na.rm=T),BRNN[L],round(time.BRNN[3]/60,digits=3))
  names(BRNN_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(BRNN_Results)
  message("")
  
  
  
  
  #cross validation for predictor_RKHS (RKHS Gaussian Kernel):
  
  message("Predicting by RKHS...")
  time.rkhs = system.time(RKHS <-runCrossVal(pheno,geno_impute,predict_RKHS,nFolds,nTimes))
  L = length(RKHS)
  message("Predict by RKHS...ended.")
  RKHS_Results <- c(mean(RKHS[1:L-1],na.rm=T),RKHS[L], round(time.rkhs[3]/60,digits=3))
  names(RKHS_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(RKHS_Results)
  message("")
  
  #cross validation for predictor_MKRKHS (MKRKHS Gaussian Kernel):
  
  message("Predicting by MKRKHS...")
  time.mkrkhs = system.time(MKRKHS <-runCrossVal(pheno,geno_impute,predict_MKRKHS,nFolds,nTimes))
  L = length(MKRKHS)
  message("Predict by MKRKHS...ended.")
  MKRKHS_Results <- c(mean(MKRKHS[1:L-1],na.rm=T),MKRKHS[L], round(time.mkrkhs[3]/60,digits=3))
  names(MKRKHS_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(MKRKHS_Results)
  message("")
  
  
  #cross validation for predictor_RF (Random Forest):
  
  message("Predicting by RF...")
  time.rf = system.time(RF <-runCrossVal(pheno,geno_impute,predict_RF,nFolds,nTimes))
  L = length(RF)
  message("Predict by RandomForest...ended.")
  RF_Results <- c(mean(RF[1:L-1],na.rm=T),RF[L],round(time.rf[3]/60,digits=3))
  names(RF_Results) = c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)")
  message(RF_Results)
  message("")
  
  
  #Table <- c(EN_Results,SVM_Results,RR_Results,LASSO_Results,rrBLUP_Results,RKHS_Results,RF_Results)
  #names(Table) <- c("Elastic-Net","SVM","RR","LASSO","rrBLUP","RKHS","RF") 
  
  Table <- cbind(BA_Results,BB_Results, BC_Results, BL_Results, BRR_Results, EN_Results,SVM_Results,BRNN_Results, RR_Results,LASSO_Results,GBLUP_Results,
                 EGBLUP_Results,RKHS_Results,MKRKHS_Results,RF_Results)
  rownames(Table) <- c("Mean Correlation","Standard Deviation","Time taken (mins)")
  colnames(Table) <- c("BayesA","BayesB","BayesC","BayesLASSO", "BayesRR", "Elastic-Net","SVM","BRNN","RR","LASSO","GBLUP","EGBLUP",
                       "RKHS","MKRKHS","RF")
  Table = t(Table) # chuyen vi
  
  return(Table)
  
}

#/////////////////////////////////////////////////////////
#THE PREDICTORS
# function to run the Genomic Prediction Models
#////////////////////////////////////////////////////////

#//////////////////////////////////////////////////////////
#   BA  Bayes A
#  uses BLLR library
#/////////////////////////////////////////////////////////////



#////////////////////////////////////////////////////////
#predict_BA uses BGLR library
########################################################################################################################

predict_BA <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  #(c)2015 V.G.TRAN & G.CHARMET
  #BayesA
  #library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) 
  {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  
  GENO = rbind(genoTrain,genoPred)
  y<- as.vector(GENO[,1])
  names(y)<- rownames(GENO)
  y[names(phenoTrain)]<- phenoTrain
  
  #whichNa<-sample(1:length(y),size=72,replace=FALSE)
  yNa<-y
  yNa[rownames(genoPred)]<-NA
  nIter=5000;
  burnIn=1000;
  thin=3;
  saveAt ='';
  S0=NULL;
  weights=NULL;
  R2=0.5;
  ETA<-list(list(X=GENO,model='BayesA'))
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(FIX)
    ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BayesA'))
    options(warn=-1)
    MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    fm=MODEL
    
  }
  
  
  else
  {
    options(warn=-1)
    #MODEL=do.call(BGLR(y=yNa,arg.BGLR))
    MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    fm=MODEL
    #MODEL1=BGLR(y=yNa,ETA=NULL,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    #MODEL2=BGLR(y=y,ETA=NULL,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  }
  
  
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  #CD = round(sqrt(1-fm$ETA1[[1]]$SD.u^2/fm$ETA1[[1]]$varU),digits=2) ; #Accuracy by coefficient of determination
  #CD = sqrt(1-fm1$SD.u^2/fm1$varU)
  CD=sqrt(1-(SDU^2/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  return(GPRED)
  
}



#////////////////////////////////////////////////////////////////////////////
# predict_BB uses BGLR library
#//////////////////////////////////////////////////////////////////////////////////////
predict_BB <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  
  library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  
  GENO = rbind(genoTrain,genoPred)
  y <- as.vector(GENO[,1])
  names(y) <- rownames(GENO)
  y[names(phenoTrain)] <- phenoTrain
  
  #whichNa<-sample(1:length(y),size=72,replace=FALSE)
  yNa <-y
  yNa[rownames(genoPred)]<-NA
  
  nIter=5000
  burnIn=1000
  thin=10
  saveAt=''
  S0=NULL
  weights=NULL
  R2=0.5
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(FIX)
    ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BayesB'))
    options(warn=-1)
    MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    fm=MODEL
    
  }
  
  else{
    ETA<-list(list(X=GENO,model='BayesB',probIn=0.05))
    options(warn=-1)
    # MODEL=do.call(BGLR,args=c(y=yNa,args.BGLR))
    MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  }
  
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  #CD = round(sqrt(1-fm$ETA1[[1]]$SD.u^2/fm$ETA1[[1]]$varU),digits=2) ; #Accuracy by coefficient of determination
  #CD = sqrt(1-fm1$SD.u^2/fm1$varU)
  CD=sqrt(1-(SDU^2/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  return(GPRED)
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_BC    uses BGLR library
#//////////////////////////////////////////////////////////////////////////////////////


predict_BC <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
  
{
  #(c)2015 V.G.TRAN & G.CHARMET
  #BayesC for Genomic Selection
  #  library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  #phenoTrain
  #phenoReal
  #genoTrain
  #genoPred
  #phenoTrain = pheno[1:250]
  #phenoReal = pheno[251:322]
  #DArT1 = MNI(DArT)
  #genoTrain = DArT1[1:250,]
  #genoPred = DArT1[251:322,]
  
  GENO <- rbind(genoTrain,genoPred)
  y <- as.vector(GENO[,1])
  names(y) <- rownames(GENO)
  y[names(phenoTrain)] <- phenoTrain
  
  #whichNa<-sample(1:length(y),size=72,replace=FALSE)
  yNa <- y
  yNa[rownames(genoPred)] <- NA
  
  nIter=5000
  burnIn=1000
  thin=3
  saveAt=''
  S0=NULL
  weights=NULL
  R2=0.5
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(FIX)
    ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BayesC'))
    options(warn=-1)
    MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    fm=MODEL
  }
  
  else
  {
    ETA<-list(list(X=GENO,model='BayesC'))
    options(warn=-1)
    #MODEL=do.call(BGLR,args=c(y=yNa,arg.BGLR))
    MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  }
  
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  #CD = round(sqrt(1-fm$ETA1[[1]]$SD.u^2/fm$ETA1[[1]]$varU),digits=2) ; #Accuracy by coefficient of determination
  #CD = sqrt(1-fm1$SD.u^2/fm1$varU)
  CD=sqrt(1-(SDU^2/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  
  return(GPRED)
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_BL Bayesian LASSO  uses BGLR library
#//////////////////////////////////////////////////////////////////////////////////////


predict_BL <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  #(c)2015 V.G.TRAN & G.CHARMET
  #BayesA
  #library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) 
  {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  
  GENO = rbind(genoTrain,genoPred)
  y<- as.vector(GENO[,1])
  names(y)<- rownames(GENO)
  y[names(phenoTrain)]<- phenoTrain
  
  #whichNa<-sample(1:length(y),size=72,replace=FALSE)
  yNa<-y
  yNa[rownames(genoPred)]<-NA
  nIter=5000;
  burnIn=1000;
  thin=3;
  saveAt='';
  S0=NULL;
  weights=NULL;
  R2=0.5;
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(FIX)
    ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BL'))
    options(warn=-1)
    MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    fm=MODEL
  }
  
  else{
    ETA<-list(list(X=GENO,model='BL'))
    
    options(warn=-1)
    #MODEL=do.call(BGLR(y=yNa,arg.BGLR))
    MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  }
  
  
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  #CD = round(sqrt(1-fm$ETA1[[1]]$SD.u^2/fm$ETA1[[1]]$varU),digits=2) ; #Accuracy by coefficient of determination
  #CD = sqrt(1-fm1$SD.u^2/fm1$varU)
  CD=sqrt(1-(SDU^2/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  
  return(GPRED)
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_BRR Bayesian rideg regression  uses BGLR library
#//////////////////////////////////////////////////////////////////////////////////////

predict_BRR <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  #(c)2015 V.G.TRAN & G.CHARMET
  #BayesA
  # library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) 
  {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  GENO = rbind(genoTrain,genoPred)
  y<- as.vector(GENO[,1])
  names(y)<- rownames(GENO)
  y[names(phenoTrain)]<- phenoTrain
  
  #whichNa<-sample(1:length(y),size=72,replace=FALSE)
  yNa<-y
  yNa[rownames(genoPred)]<-NA
  nIter=8000;
  burnIn=2
  
  
  000;
  thin=3;
  saveAt='';
  S0=NULL;
  weights=NULL;
  R2=0.5;
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(FIX)
    ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BRR'))
    options(warn=-1)
    MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    fm=MODEL
  }
  
  else
  {
    ETA<-list(list(X=GENO,model='BRR'))
    
    options(warn=-1)
    #MODEL=do.call(BGLR(y=yNa,arg.BGLR))
    MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  }
  
  
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  #CD = round(sqrt(1-fm$ETA1[[1]]$SD.u^2/fm$ETA1[[1]]$varU),digits=2) ; #Accuracy by coefficient of determination
  #CD = sqrt(1-fm1$SD.u^2/fm1$varU)
  CD=sqrt(1-(SDU^2/VARU))
  CD=1-(SDU^2/VARU)
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  
  return(GPRED)
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_ElasticNet  uses glmnet library
#//////////////////////////////////////////////////////////////////////////////////////


# predict_ElasticNet <- function(phenoTrain, genoTrain, genoPred,arg.cv.glmnet=list(family="gaussian",alpha=0.5)){
predict_ElasticNet <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  
  #(c)2013 V.G. TRAN
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  #  library(glmnet)  # LASSO & Elastic-Net generalized linear models
  
  #do cross-validation to get the optimal value of lambda:
  
  # cv.fit <- do.call(glmnet::cv.glmnet,args=c(genoTrain,phenoTrain,arg.cv.glmnet)) #ElasticNet penalty with top results
  cv.fit <- cv.glmnet(genoTrain,phenoTrain,family="gaussian",alpha=0.5) #ElasticNet penalty with top results
  
  #alpha=1: lasso penalty, alpha=0: ridge penalty
  lambda_min <- cv.fit$lambda.min # making the best prediction
  #lambda.1se <- cv.fit$lambda.1se
  
  gpred <- predict(cv.fit,newx=genoPred,s=c(lambda_min))
  #GPRED <- round(gpred,digits=3)
  gpredSD=gpred
  gpredSD=NA
  CD <- gpred
  CD <- NA
  GPRED=cbind(gpred,gpredSD,CD)
  GPRED=round(GPRED,digits=3)
  return(GPRED)
}

#////////////////////////////////////////////////////////////////////////////
# predict_BRNN Bayesian Regularization Neural Network  uses BRNN library
#//////////////////////////////////////////////////////////////////////////////////////


predict_BRNN<- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  
  #(c)2013 V.G. TRAN
  if (!requireNamespace("brnn", quietly = TRUE)) {
    stop("brnn needed for this function to work. Please install it.",
         call. = FALSE)
  }
  library(brnn)
  #  library(glmnet)  # LASSO & Elastic-Net generalized linear models
  
  #do cross-validation to get the optimal value of lambda:
  
  # cv.fit <- do.call(glmnet::cv.glmnet,args=c(genoTrain,phenoTrain,arg.cv.glmnet)) #ElasticNet penalty with top results
  cv.fit <- brnn(genoTrain,phenoTrain,neurons=2,epochs=10, verbose=T) # neural networks with 2 neurons and 50 iterations
  
  
  gpred <- predict(cv.fit,newdata=genoPred)
  #GPRED <- round(gpred,digits=3)
  gpredSD=gpred
  gpredSD=NA
  CD <- gpred
  CD <- NA
  GPRED=cbind(gpred,gpredSD,CD)
  GPRED=round(GPRED,digits=3)
  rownames(GPRED)=rownames(genoPred)
  return(GPRED)
}

#////////////////////////////////////////////////////////////////////////////
# predict_GBLUP  uses  rrBLUP 
#//////////////////////////////////////////////////////////////////////////////////////


predict_GBLUP <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred) 
  # NB based on kin.blu (rrBLUP), faster than bayesian, independant from marker number
{
  
  
  #BLUP #(c)2015 V.G. TRAN & G. CHARMET
  #  library(BGLR)
  #       if (!requireNamespace("rrBLUP", quietly = TRUE)) {
  #   stop("rrBLUP needed for this function to work. Please install it.",
  #    call. = FALSE)
  #  }
  
  
  GENO <- rbind(genoTrain,genoPred)
  
  
  y <- as.vector(GENO[,1])
  names(y) <- rownames(GENO)
  y[names(phenoTrain)] <- phenoTrain
  yNa <- y
  yNa[rownames(genoPred)] <- NA
  A <-A.mat(GENO) 
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(MNI(FIX))
    #ixed=colnames(FIX)
    dataF=data.frame(genoID=names(y),yield=yNa,FIX)
    Fixed=colnames(dataF[,-c(1,2)])
    
    MODEL=kin.blup(dataF,geno="genoID",pheno="yield", GAUSS=FALSE,K=A, fixed=Fixed,
                   PEV=TRUE,n.core=1,theta.seq=NULL)
  }
  else
  {
    
    dataF=data.frame(genoID=names(y),yield=yNa)
    #Fixed=colnames(dataF[,-c(1,2)])
    
    MODEL=kin.blup(data=dataF,geno="genoID",pheno="yield", GAUSS=FALSE, K=A, 
                   PEV=TRUE,n.core=1,theta.seq=NULL)
  }
  
  
  
  options(warn=-1)
  # PREDICT <- do.call(BGLR::BGLR,args=c(y=yNa,ETA=ETA,arg.BGLR))
  
  
  gpred = MODEL$pred
  gpredSD=sqrt(MODEL$PEV)
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=MODEL$Vg   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$PEV
  
  
  #CD: coefficient of determination:
  
  
  CD=sqrt(1-(SDU/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  
  return(GPRED)
  
}

predict_GBLUPB <- function(phenoTrain, genoTrain, genoPred,arg.kinship.BLUP=list(K.method="GAUSS")) 
  # OLD function based on BGLR: too long with bayesian solution
{
  
  
  #BLUP #(c)2015 V.G. TRAN & G. CHARMET
  #  library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  GENO <- rbind(genoTrain,genoPred)
  y <- as.vector(GENO[,1])
  names(y) <- rownames(GENO)
  y[names(phenoTrain)] <- phenoTrain
  yNa <- y
  yNa[rownames(genoPred)] <- NA
  
  A <-A.mat(GENO) 
  
  #2# Computing D and then K
  
  
  #3# Kernel Averaging using BGLR
  # ETA <- list(list(K=exp(-h[1]*D),model='RKHS'),
  #          list(K=exp(-h[2]*D),model='RKHS'),
  #         list(K=exp(-h[3]*D),model='RKHS'))
  ETA= list(list(K=A, model='RKHS'))
  # options(warn=-1)
  # PREDICT <- do.call(BGLR::BGLR,args=c(y=yNa,ETA=ETA,arg.BGLR))
  MODEL <- BGLR(y=yNa,ETA=ETA,nIter=5000, burnIn=1000)
  
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  
  CD=sqrt(1-(SDU^2/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  
  return(GPRED)
  
}





#////////////////////////////////////////////////////////////////////////////
# predict_MKRKHS  Multiple Kernel Reproductive Kernel Hilbert Space uses BGLR library 
#//////////////////////////////////////////////////////////////////////////////////////


# predict_MKRKHS <- function(phenoTrain, genoTrain, genoPred,arg.BGLR=list(nIter=5000, burnIn=1000)) {
predict_MKRKHS <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred) 
{
  #Multi-Kernel RKHS
  #(c)2015 V.G. TRAN & G. CHARMET
  #  library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  GENO <- rbind(genoTrain,genoPred)
  y <- as.vector(GENO[,1])
  names(y) <- rownames(GENO)
  y[names(phenoTrain)] <- phenoTrain
  yNa <- y
  yNa[rownames(genoPred)] <- NA
  
  
  #A = A.mat(GENO);
  X = GENO
  p = ncol(X)
  
  #2# Computing D and then K
  X <- scale(X,center=TRUE,scale=TRUE)
  D <- (as.matrix(dist(X,method='euclidean'))^2)/p
  # h<-0.5*c(1/5,1,5)
  
  h <- sqrt(c(0.2,0.5,0.8))
  
  #3# Kernel Averaging using BGLR
  ETA <- list(list(K=exp(-h[1]*D),model='RKHS'),
              list(K=exp(-h[2]*D),model='RKHS'),
              list(K=exp(-h[3]*D),model='RKHS'))
  
  options(warn=-1)
  # PREDICT <- do.call(BGLR::BGLR,args=c(y=yNa,ETA=ETA,arg.BGLR))
  MODEL <- BGLR(y=yNa,ETA=ETA,nIter=8000, burnIn=2000)
  
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  #CD = round(sqrt(1-fm$ETA1[[1]]$SD.u^2/fm$ETA1[[1]]$varU),digits=2) ; #Accuracy by coefficient of determination
  #CD = sqrt(1-fm1$SD.u^2/fm1$varU)
  CD=sqrt(1-(SDU^2/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  
  return(GPRED)
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_RKHS Reproductive Kernel Hilbert Space uses BGLR library 
#//////////////////////////////////////////////////////////////////////////////////////


predict_RKHS <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred) 
{
  # RKHS
  #(c)2015 V.G. TRAN & G. CHARMET
  #  library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  GENO <- rbind(genoTrain,genoPred)
  y <- as.vector(GENO[,1])
  names(y) <- rownames(GENO)
  y[names(phenoTrain)] <- phenoTrain
  yNa <- y
  yNa[rownames(genoPred)] <- NA
  
  X = GENO
  p = ncol(X)
  
  #2# Computing D and then K
  X <- scale(X,center=TRUE,scale=TRUE)
  D <- (as.matrix(dist(X,method='euclidean'))^2)/p
  h <- 0.5
  K=exp(-h*D)
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(FIX)
    ETA2<-list(list(X=FIX,model="FIXED"),list(K=K, model='RKHS'))
    options(warn=-1)
    MODEL=BGLR(y=yNa,ETA=ETA2,nIter=5000, burnIn=1000)
    fm=MODEL
  }
  else
  {
    
    ETA= list(list(K=K, model='RKHS'))
    options(warn=-1)
    # PREDICT <- do.call(BGLR::BGLR,args=c(y=yNa,ETA=ETA,arg.BGLR))
    MODEL <- BGLR(y=yNa,ETA=ETA,nIter=5000, burnIn=1000)
  }
  gpred = MODEL$yHat
  gpredSD=MODEL$SD.yHat
  GPRED=cbind(gpred,gpredSD)
  
  
  VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
  SDU=MODEL$SD.yHat
  
  
  #CD: coefficient of determination:
  
  #CD = round(sqrt(1-fm$ETA1[[1]]$SD.u^2/fm$ETA1[[1]]$varU),digits=2) ; #Accuracy by coefficient of determination
  #CD = sqrt(1-fm1$SD.u^2/fm1$varU)
  CD=sqrt(1-(SDU^2/VARU))
  GPRED=cbind(GPRED,CD)
  GPRED=round(GPRED,digits=3)
  GPRED=GPRED[rownames(genoPred),]
  
  return(GPRED)
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_RR Random Forest regression uses randomForest library 
#//////////////////////////////////////////////////////////////////////////////////////


predict_RF <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred) 
{
  #(c)2013 V.G.TRAN
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("randomForest needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  #  library(randomForest)
  #return(randomForest(genoTrain, phenoTrain, xtest=genoPred)$test$predicted)
  #gpred <- do.call(randomForest::randomForest,args=c(genoTrain, phenoTrain, xtest=genoPred,arg.randomForest))$test$predicted
  gpred = randomForest(genoTrain, phenoTrain, xtest=genoPred)$test$predicted;
  gpredSD=gpred
  gpredSD=NA
  CD <- gpred
  CD <- NA
  GPRED=cbind(gpred,gpredSD,CD)
  GPRED=round(GPRED,digits=3)
  return(GPRED)
  
}


#//////////////////////////////////////////////////////////////////////////////////
#predict_EGBLUP Epistatic GBLUP uses BGLR library
#//////////////////////////////////////////////////////////////////////////////////

predict_EGBLUP <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  #(c)2015 G.CHARMET & V.G. TRAN
  #EG-BLUP for Genomic Selection
  #library(BGLR)
  
  GENO = rbind(genoTrain,genoPred)
  A_GENO <- AM(GENO)
  y<- as.vector(GENO[,1])
  names(y)<- rownames(GENO)
  y[names(phenoTrain)]<- phenoTrain
  
  #whichNa<-sample(1:length(y),size=72,replace=FALSE)
  yNa<-y
  yNa[rownames(genoPred)]<-NA
  
  nIter=5000;
  burnIn=1000;
  thin=3;
  saveAt='';
  S0=NULL;
  weights=NULL;
  R2=0.5;
  
  if(FixedTrain!="NULL")
  {
    FIX <- rbind(FixedTrain,FixedPred)
    FIX=round(FIX)
    ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BRR'),list(K=A_GENO*A_GENO,model='RKHS'))
    options(warn=-1)
    MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
    fm=MODEL
  }
  
  
  else{
    ETA<-list(list(X=GENO,model='BRR'),list(K=A_GENO*A_GENO,model='RKHS'))
    options(warn=-1);
    MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  }
  
  gpred = MODEL$yHat
  #gpred=predict(MODEL)
  gpredSD=MODEL$SD.yHat
  VARU=var(gpred)
  SDU=gpredSD
  CD = round(sqrt(1-(SDU^2/VARU)),digits=3) ; #Accuracy by coefficient of determination
  GPRED=cbind(gpred,gpredSD,CD)
  GPRED=GPRED[rownames(genoPred),]
  
  
  return(GPRED)
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_RR Ridge regression uses glmnet library 
#//////////////////////////////////////////////////////////////////////////////////////


predict_RR <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred) 
{
  
  #(c)2013 V.G.TRAN
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet needed for this function to work. Please install it.",
         call. = FALSE)
  }  
  
  #  library(glmnet) 
  
  #do cross-validation to get the optimal value of lambda:
  # cv.fit <- do.call(glmnet::cv.glmnet,args=c(genoTrain,phenoTrain,arg.cv.glmnet)) #Ridge penalty with top results
  cv.fit <- cv.glmnet(genoTrain,phenoTrain,alpha=0) #Ridge penalty with top results
  
  lambda_min <- cv.fit$lambda.min # making the best prediction
  #return(predict(cv.fit,newx=genoPred,s=c(lambda_min)))
  
  gpred <- predict(cv.fit,newx=genoPred,s=c(lambda_min))
  gpredSD=gpred
  gpredSD=NA
  CD <- gpred
  CD <- NA
  GPRED=cbind(gpred,gpredSD,CD)
  GPRED=round(GPRED,digits=3)
  return(GPRED)
  
  
}

#////////////////////////////////////////////////////////////////////////////
# predict_Lasso LASSO uses glmnet library 
#//////////////////////////////////////////////////////////////////////////////////////

predict_Lasso <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred)
{
  #(c)2013 V.G.TRAN
  
  library(glmnet) 
  
  #do cross-validation to get the optimal value of lambda:
  
  cv.fit <- cv.glmnet(genoTrain,phenoTrain,family="gaussian",alpha=1) #Lasso penalty with top results
  
  lambda_min <- cv.fit$lambda.min # making the best prediction
  
  gpred = predict(cv.fit,newx=genoPred,s=c(lambda_min))
  gpredSD=gpred
  gpredSD=NA
  CD <- gpred
  CD <- NA
  GPRED=cbind(gpred,gpredSD,CD)
  GPRED=round(GPRED,digits=3)
  return(GPRED) 
}

#////////////////////////////////////////////////////////////////////////////
# predict_SVM Support vector machine uses e1071 library 
#//////////////////////////////////////////////////////////////////////////////////////


predict_SVM <- function(phenoTrain, genoTrain, FixedTrain, genoPred, FixedPred) 
{
  
  #(c)2013 vangiang.tran@gmail.com
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("e1071 needed for this function to work. Please install it.",
         call. = FALSE)
  }  
  
  #  library(e1071)
  #model <- do.call(e1071::svm,args=c(genoTrain,phenoTrain,arg.svm))
  # model <- svm(genoTrain,phenoTrain,method="C-classification",kernel="radial",cost=10,gamma=0.001)
  model <- svm(genoTrain,phenoTrain,method="nu-regression",kernel="radial",cost=10,gamma=0.001)
  gpred <- predict(model,genoPred)
  gpredSD=gpred
  gpredSD=NA
  CD <- gpred
  CD <- NA
  GPRED=cbind(gpred,gpredSD,CD)
  GPRED=round(GPRED,digits=3)
  return(GPRED) 
}


#END OF PREDICTOR FUNCTIONS
#//////////////////////////////////////////////////////////////


transfer <- function(data,time.cal,ncol_geno_shrink) {
  #(c)2014 V.G. TRAN
  SUMMA <- data[[1]]
  L <- length(SUMMA)
  Summary <- c(mean(SUMMA[1:L-1],na.rm=T),SUMMA[L],round(time.cal[3]/60,digits=2),ncol_geno_shrink)
  names(Summary) <- c("Mean Correlation CVs","Standard Deviation CVs","Time taken (mins)","Number of markers")
  print(Summary)
  
  
  BV_TABLE <- data[[3]]
  CV_all_times <- SUMMA[1:L-1]
  names(CV_all_times) <- paste("time",1:length(CV_all_times))
  SD_all_times <- SUMMA[L]
  names(SD_all_times) <-"sd"
  
  MSEP <- data[[2]]
  MSEP_all_times <- MSEP[1:L-1]
  names(MSEP_all_times) <- paste("time",1:length(MSEP_all_times))
  SD_MSEP <- MSEP[L]
  names(SD_MSEP) <-"sd squarred MSEP)"
  
  Results <- list(Summary,CV_all_times,SD_all_times,MSEP_all_times,SD_MSEP,BV_TABLE)
  
  names(Results) <- c("summary","cv","SD_cv","MSEP","SD_MSEP","bv_table")
  message("")
  message("Processing...ended! Results in: object$summary, object$cv, object$sd, object$bv")
  return(Results)
}

#waitingbar <- function(x)



waitingbar <- function(x = sort(runif(20)), ...) {
  #(c)2013 V.G. TRAN
  
  pb <- txtProgressBar(min=0,max=1,initial=0,char="*",width=40,style=3)
  for(i in c(0, x, 1)) {Sys.sleep(0.1)
    setTxtProgressBar(pb, i)}
  Sys.sleep(0.1)
  close(pb)
}

#To use:
#waiting.bar()
#waiting.bar(runif(10))
#waiting.bar(style = 3)

#////////////////////////////////////////////////////////////////////////////
#QTL SIMULATION:
#//////////////////////////////////////////////////////////////////////////////////////

qtlSIM <- function(geno,NQTL=100,h2=0.3) {
  #(c)2014 G. CHARMET & V.G. TRAN
  
  
  
  RealizedH2 <- numeric()  # realized trait heritability: to check if it fits expected value
  
  distriQTLh2 <- numeric()# realized QTL h?, to plot histogrammes....
  
  
  S <- sample(ncol(geno),NQTL)
  S <- sort.list(S)
  
  QTL <- geno[,S]
  X <- geno[,-S]
  
  # simulation of effects as a single gaussian distribution
  Effects <- rnorm(NQTL)    # Gaussian distribution
  
  QTL <- t(t(QTL)*Effects)
  
  TBV <- apply(QTL,1,sum)
  varG <- var(TBV)
  varE <- varG*((1-h2)/h2)
  
  noise <- rnorm(length(TBV),0,sqrt(varE))
  QT <- TBV+noise
  
  h2QTL <- Effects^2/var(QT)
  
  
  result <- list(newSNP=X,pheno=QT, TBV=TBV, Effects=Effects, h2QTL=h2QTL )
  
  result
  
}

#/////////////////////////////////////////////
bwgs.nacount <- function(geno) {
  #Calculate the percentage of missing elements in SNP or DArT, SNP GBS:
  #(c)06/02/2014 V.G. Tran
  
  num_element = table(geno)
  num_of_binary = length(num_element) # = 2 for DArT, = 3 for SNP, = 5 for SNP issu GBS
  num_element_total = dim(geno)[1]*dim(geno)[2]
  
  if(num_of_binary==2){
    
    num_na = num_element_total-(num_element[1]+num_element[2])
    
  } else {
    
    if(num_of_binary==3){
      
      num_na = num_element_total-(num_element[1]+num_element[2]+num_element[3])
      
    } else {
      
      if(num_of_binary==4){
        num_na = num_element_total-(num_element[1]+num_element[2]+num_element[3]+num_element[4])
      } else {
        
        if(num_of_binary==5){
          num_na = num_element_total-(num_element[1]+num_element[2]+num_element[3]+num_element[4]+num_element[5])
        } else {
          
          if(num_of_binary==6){
            num_na = num_element_total-(num_element[1]+num_element[2]+num_element[3]+num_element[4]+num_element[5]+num_element[6])
          } 
          
          else {
            stop("Error in The genotype data (or non numerized).")
          } # 
        }
      }
    }
  }
  
  na_pourcentage = num_na*100/num_element_total # in percentage
  
  na_pourcentage = as.vector(na_pourcentage)
  
  return(na_pourcentage)
  
}

#////////////////////////////////////////////////////////////
#COMMON FUNCTION
bwgs.common <- function(a1,a2) {
  #COMMON FUNCTION
  commun = a1[which(names(a1)%in%names(a2))]
  return(commun)
  
}

bwgs.bestpred <- function(ans,n) {
  #(c)2014 G. CHARMET & V.G. TRAN
  message(c("The ",as.character(n)," best prediction[s]:"))
  message("")
  #Find the n best predictions:
  ans.BEST = ans[,sort.list(-ans)]
  ans.bestn = ans.BEST[1:n]
  #message(c("The ",n," best prediction are: ",ans.bestn))
  return(ans.bestn)
}



#///////////////////////////////////////////////////////
#MNI Function replaces missing value by average allele frequency
#////////////////////////////////////////////////////////
MNI <- function(x) 
{
  NAreplace= function (y) 
  {
    if( length(na.omit(y))!=0){
      y[is.na(y)] = mean(y, na.rm = TRUE)
    }
    return(y)
  }
  imp= apply(x, 2, NAreplace)
  return(imp)
} 

#-------------------------------------------------
#Function for optimizing the calibration set
# by selecting a subset in the reference (training population)
# which maximizes the CDmean criteria (Rincent et al 2012)
#-------------------------------------------------

# Main function
# apply CDmean estimation either on rendom or optimized samples
# Infiles required are

# matA1 = the additive relationships matrix, e.g. estimated from marker data using pedigree or A.mat
# pheno = the vector of phenotypes for the training set corresponding to matA1
# h2 the heritability of the pheno trait
# Nindrep = the number of genotypes in the subsample
# otptions are:
# Method = "Boot to indicate whether random sampling (Bootstraps) should be done
#   or "Opti" to indicate that the optimization

optiTRAIN=function(geno,NSample=100,Nopti=3000)  #  supress pheno, useless
{
  # Preliminary computations to be run once
  matA1=A.mat(geno)
  matA1=as.matrix(matA1)
  
  Nind=nrow(matA1) # total number of individuals  
  #varP=var(pheno) # Pheno is a vector of phenotypes
  #h2=0.5
  
  #varG=h2*varP
  #varE=(1-h2)/h2*varG
  #lambda=varE/varG # lambda is needed to estimate the CDmean
  lambda=1
  invA1=solve(matA1)
  
  
  ##############################
  # Optimization algorithm
  ##############################
  
  
  #Design matrices
  Ident<-diag(NSample)
  X<-rep(1,NSample)
  M<-Ident- (X%*%solve(t(X)%*%X) %*% t(X) )
  
  Sample1<-sample(Nind,NSample) #Calibration set initialization
  SaveSample=Sample1
  NotSampled1<-seq(1:Nind)
  NotSampled<-NotSampled1[-Sample1] # Initial validation set
  
  Z=matrix(0,NSample,Nind)
  for (i in 1:length(Sample1)) { Z[i,Sample1[i]]=1 } 
  
  T<-contrasteNonPheno(NotSampled,Nind,NSample)   # T matrice des contrastes
  
  # Calculate of CDmean of the initial set
  matCD<-(t(T)%*%(matA1-lambda*solve(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%matA1%*%T)
  CD=diag(matCD)
  CDmeanSave=mean(CD)
  
  CDmeanMax1=rep(NA,Nopti)
  
  # Exchange algorithm (maximize CDmean)
  cpt2=1
  cpt=0
  while (cpt2<Nopti) {  # Make sure that Nopti is enough in your case (that you reached a plateau), for this look at CDmeanMax1.
    NotSampled=NotSampled1[-Sample1] 
    cpt2=cpt2+1
    # Remove one individual (randomly choosen) from the sample :
    Sample2=sample(Sample1,1)
    # Select one individual (randomly choosen) from the individuals that are not in the Calibration set :
    Sample3=sample(NotSampled,1)
    # New calibration set :
    Sample4=c(Sample3,Sample1[Sample1!=Sample2])
    # Calculate the mean CD of the new calibration set :
    Z=matrix(0,NSample,Nind)
    for (i in 1:length(Sample4)) { Z[i,Sample4[i]]=1 } 
    NotSampled=NotSampled1[-Sample4] 
    T<-contrasteNonPheno(NotSampled,Nind,NSample)
    
    matCD<-(t(T)%*%(matA1-lambda*solve(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%matA1%*%T)
    CD=diag(matCD)
    
    if (mean(CD)>CDmeanSave ) { Sample1=Sample4 # Accept the new Calibration set if CDmean is increased, reject otherwise.
    CDmeanSave=mean(CD)  
    cpt=0 } else { cpt=cpt+1 
    }
    CDmeanMax1[cpt2-1]=CDmeanSave
  }  #Fin du while
  
  sampleOptimiz=Sample1 # SampleOptimiz is the optimized calibration set
  sampleOptimiz=Sample1 # SampleOptimiz is the optimized calibration set
  sortOptimiz=sampleOptimiz[sort.list(sampleOptimiz)]
  genoOptimiz=geno[sortOptimiz,]
  result=list(CDmax=CDmeanSave,sampleOPTI=sampleOptimiz,genoOptimiz=genoOptimiz)
  
  
}
# This function creates the matrix of contrast between each of the individual not in the calibration set and the mean of the population
contrasteNonPheno=function(NotSampled,Nind,NSample)
{
  mat=matrix(-1/Nind,Nind,Nind-NSample)
  for (i in 1:ncol(mat)) {
    mat[NotSampled[i],i]=1-1/Nind
  }
  return(mat)
}

