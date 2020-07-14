#' Computes the GEBV prediction for the target population with only genotypic Data using the options for model selection.
#'
#' @param geno_train  Matrix (n x m) of genotypes for the training population: n lines with m markers. Genotypes should be coded as {-1, 0, 1, NA}. Missing data are allowed and coded as NA.
#' @param pheno_train Vector (n x 1) of phenotype for the training phenotypes. This vector should have no missing values. Otherwise, missing values (NA) will be omitted in both pheno_train and geno_train.
#' @param geno_target Matrix (z x m) of genotypes for the target population: z lines with the same m markers as in geno_train. Genotypes should be coded as {-1, 0, 1, NA}. Missing data are allowed and coded as NA. Other arguments are identical to those of bwgs.cv, except pop_reduct_method, nTimes and nFolds, since the prediction is run only once, using the whole training population for model estimation, then applied to the target population.
#' @param FIXED_train A matrix of fixed effect for training, to be used with some methods such as those included in BGLR, MUST have same rownames as geno and coded(-1 0 1)
#' @param FIXED_target A matrix of fixed effect for targeting, to be used with some methods such as those included in BGLR, MUST have same rownames as geno and coded(-1 0 1)
#' @param MAXNA The maximum proportion of missing value which is admitted for filtering marker columns in geno. Default value is 0.2
#' @param MAF The minimum allele frequency for filtering marker colums in geno; default value is 0.05
#' @param geno.reduct.method Allows sampling a subset of markers for speeding up computing time and/or avoid introducing more noise than informative markers. Options are:
#' \itemize{
#'     \item{RMR: Random sampling (without replacement) of a subset of markers. To be used with the parameter “reduct.marker.size”.}
#'     \item{LD (with r2 and MAP): enables “pruning” of markers which are in LD > r2. Only the marker with the least missing values is kept for each pair in LD>r2. To allow faster computation, r2 is estimated chromosome by chromosome, so a MAP file is required with information of marker assignation to chromosomes. The MAP file should contain at least three columns: marker_name, chromosome_name and distance_from_origin (either genetic of physical distance, only used for sorting markers, LD being re-estimated from marker Data).}
#'     \item{ANO (with pval): one-way ANOVA are carried out with R function lm on trait “pheno”. Every markers are tested one at a time, and only markers with pvalue<pval are kept for GEBV prediction.}
#'     \item{ANO+LD (with pval and r2, MAP is facultative): combines a first step of marker selection with ANO, then a second step of pruning using LD option.}
#' }
#' @param reduct.size Specifies the number of markers for the genotypic reduction using RMR (reduct.size < m).
#' @param r2 Coefficient of linkage disequilibrium (LD). Setting 0<r2<1 if the genotypic reduction method is in {LD or ANO+LD }.
#' @param pval p value for ANO method, 0 < pval < 1.
#' @param MAP A file with markers in rows ane at least ONE columns with colnames= "chrom". Used for computing r2 within linkage groups.
#' @param geno.impute.method Allow missing marker data imputation using the two methods proposed in function A.mat of package rrBLUP, namely:
#' \itemize{
#'     \item{MNI: missing data are replaced by the mean allele frequency of the marker (column in geno)}
#'     \item{EMI: missing data are replaced using an expectation-maximization methods described in function A.mat (Endelman & Janninck 2012).}
#' }
#' 
#' Default value is NULL.
#' 
#' Note that these imputation methods are only suited when there are a few missing value, typically in marker data from SNP chips of KasPAR. They are NOT suited for imputing marker data from low density to high density designs, and when there are MANY missing Data as typically provided by GBS. More sophisticated software (e.g. Beagles, Browning & Browning 2016) should be used before BWGS.
#' @param predict.method The options for genomic breeding value prediction methods. The available options are:
#' \itemize{
#'    \item{GBLUP: performs G-BLUP using a marker-based relationship matrix, implemented through BGLR R-library. Equivalent to ridge regression (RR-BLUP) of marker effects.}
#'    \item{EGBLUP: performs EG-BLUP, i.e. BLUP using a "squared" relationship matrix to model epistatic 2x2 interactions, as described by Jiang & Reif (2015), using BGLR library}
#'    \item{RR: ridge regression, using package glmnet. In theory, strictly equivalent to gblup.}
#'    \item{LASSO: Least Absolute Shrinkage and Selection Operator is another penalized regression methods which yield more shrinked estimates than RR. Run by glmnet library.}
#'    \item{EN: Elastic Net (Zou and  Hastie, 2005), which is a weighted combination of RR and LASSO, using glmnet library}
#'  }
#'  Several Bayesian methods, using the BGLR library:
#'  \itemize{
#'      \item{BRR: Bayesian ridge regression: same as rr-blup, but bayesian resolution. Induces homogeneous shrinkage of all markers effects towards zero with Gaussian distribution (de los Campos et al, 2013)}
#'      \item{BL: Bayesian LASSO: uses an exponential prior on marker variances priors, leading to double exponential distribution of marker effects (Park & Casella 2008)}
#'      \item{BA: Bayes A uses a scaled-t prior distribution of marker effects. (Meuwissen et al 2001).}
#'      \item{BB: Bayes B, uses a mixture of distribution with a point mass at zero and with a slab of non-zero marker effects with a scaled-t distribution (Habier et al 2011).}
#'      \item{BC: Bayes C same as Bayes B with a slab with Gaussian distribution.}
#' }
#' A more detailed description of these methods can be found in Perez & de los Campos 2014 (http://genomics.cimmyt.org/BGLR-extdoc.pdf).
#' Three semi-parametric methods:
#' \itemize{
#'     \item{RKHS: reproductive kernel Hilbert space and multiple kernel MRKHS, using BGLR (Gianola and van Kaam 2008).  Based on genetic distance and a kernel function to regulate the distribution of marker effects. This methods is claimed to be effective for detecting non additive effects.}
#'     \item{RF: Random forest regression, using randomForest library (Breiman, 2001, Breiman and  Cutler 2013). This methods uses regression models on tree nodes which are rooted in bootstrapping data. Supposed to be able to capture interactions between markers}
#'     \item{SVM: support vector machine, run by e1071 library. For details, see Chang, Chih-Chung and Lin, Chih-Jen: LIBSVM: a library for Support Vector Machines http://www.csie.ntu.edu.tw/~cjlin/libsvm}
#'     \item{BRNN: Bayesian Regularization for feed-forward Neural Network, with the R-package BRNN (Gianola et al 2011). To  keep computing time in reasonable limits, the parameters for the brnn function are neurons=2 and epochs = 20.}
#' }
#'
#' @return
#' The object bwgs.predict returns Matrix of dimension nx3. Columns are:
#' \itemize{
#'     \item{Predict BV: the nx1 vector of GEBVs for the validation set (rows of geno_valid)}
#'     \item{gpredSD: Standart deviation of estimated GEBV}
#'     \item{CD: coefficient of determination for each GEBV, estimated as sqrt ((1-stdev(GEBVi))^2/2g)}
#' }
#' Note that gpredSD and CD are only available for methods using the BGLR library, namely GBLUP, EGBLUP, BA,BB,BC,BL,RKHS and MKRKHS. 
#' These two columns contain NA for methods RF, RR, LASSO, EN and SVM.
#' @examples 
#' \dontrun{
#' data(inra)
#' testPREDICT_GBLUP <- bwgs.predict(geno_train = TRAIN47K,
#'      pheno_train = YieldBLUE,
#'      geno_target = TARGET47K,
#'      MAXNA = 0.2, 
#'      MAF = 0.05,
#'      geno.reduct.method = "NULL",
#'      reduct.size = "NULL",
#'      r2 = "NULL",
#'      pval = "NULL",
#'      MAP = "NULL",
#'      geno.impute.method = "MNI",
#'      predict.method = "GBLUP")


#' }
#' @export
bwgs.predict <- function(geno_train,pheno_train,geno_target,FIXED_train="NULL",FIXED_target="NULL",MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",r2="NULL",pval="NULL",
                         MAP="NULL",geno.impute.method="NULL",predict.method="GBLUP") 
{
  #(c)2015 louis.gautier.tran@gmail.com & gilles.charmet@clermont.inra.fr
  message("2017 BWGS  - Version 1.6.0 Release date: 31/02/2017")
  #message("2015 Gilles Charmet & Louis Gautier Tran)
  
  # upload the necessary libraries

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
      time.chrld = system.time(geno_shrink <- CHROMLD(geno_shrinkANO,R2seuil = r2,MAP2))
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
      time.chrld = system.time(geno_shrink <- CHROMLD(geno_impute,R2seuil=r2,MAP))
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


