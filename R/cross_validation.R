
#Title  bwgs.cv
#' @title Genomic Prediction with cross validation
#' @description The bwgs.cv function carries out cross-validation using genotypic and
#' phenotypic data from a reference population, with options for genotypic
#' matrix processing and genomic breeding value estimation.
#' @param geno Matrix (n x m) of genotypes for the training population: n lines with m markers. Genotypes should be coded {-1, 0, 1}. Missing data are allowed and coded as NA.
#' @param pheno Vector (n x 1) of "phenotypes", i.e. observations or pre-processed, corrected values. This vector should have no missing values, otherwise missing values (NA) will be omitted in both pheno and geno. In a first step, bwgs.cv checks  whether rownames(geno) match with  names(pheno). If not the case,  the common elements (intersect)  are selected in both geno and  pheno for further analyses. If a  MAP file is provided, the  selected set of markers are also  sorted out in MAP.
#' @param FIXED A matrix of fixed effect, to be used with some methods such as those included in BGLR, MUST have same rownames as geno and coded(-1 0 1)
#' @param MAXNA The maximum proportion of missing  value which is admitted for filtering marker columns in geno. Default value is 0.2
#' @param MAF The minimum allele frequency for filtering marker colums in geno;default value is 0.05
#' @param pop.reduct.method Method for reducing the size of the training population. Can be used for teaching purposes, no real interest in real life if the entire population is already genotyped and phenotyped. 
#' Default value is NULL (all training set used).
#' Proposed methods are:
#' \itemize{
#'    \item{RANDOM: a subset of sample.pop.size is randomly selected for training the model, and the unselected part of the population is used for validation. The process is repeated nFolds *  nTimes to have the same number of replicates than with cross-validation.}
#'    \item{OPTI: the optimization algorithm based on CDmean (Rincent et al 2012) to select a subset which maximizes average CD (coefficient of determination) in the validation set. Since the process is long and has some  stochastic component, it is  repeated only nTimes.}
#' }
#' @param sample.pop.size The size of the subset of individuals in the training set (both geno and pheno) selected by pop.reduct.method if not NULL.
#' @param geno.reduct.method Allows sampling a subset of markers for speeding up computing time and/or avoid introducing more noise than informative markers. Options are:
#' \itemize{
#'     \item{RMR: Random sampling (without replacement) of a subset of markers. To be used with the parameter “reduct.marker.size”.}
#'     \item{LD (with r2 and MAP): enables “pruning” of markers which are in LD > r2. Only the marker with the least missing values is kept for each pair in LD>r2. To allow faster computation, r2 is estimated chromosome by chromosome, so a MAP file is required with information of marker assignation to chromosomes. The MAP file should contain at least three columns: marker_name, chromosome_name and distance_from_origin (either genetic of physical distance, only used for sorting markers, LD being re-estimated from marker Data).}
#'     \item{ANO (with pval): one-way ANOVA are carried out with R function lm on trait “pheno”. Every markers are tested one at a time, and only markers with pvalue<pval are kept for GEBV prediction.}
#'     \item{ANO+LD (with pval and r2, MAP is facultative): combines a first step of marker selection with ANO, then a second step of pruning using LD option.}
#' }
#' @param reduct.marker.size Specifies the number of markers for the genotypic reduction using RMR (reduct.size < m).
#' @param pval p value for ANO method, 0 < pval < 1.
#' @param r2 Coefficient of linkage disequilibrium (LD). Setting 0<r2<1 if the genotypic reduction method is in {LD or ANO+LD }.
#' @param MAP A matrix with markers in rows and at least ONE columns with colnames= "chrom". Used for computing r2 within linkage groups.
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
#' @param nFolds Number of folds for the cross-validation. Smallest value recommended is nFolds = 3.
#' @param nTimes Number of independent replicates for the cross-validation. Smallest value recommended is nTimes = 3.
#' @return 
#' The class bwgs.cv returns a list containing:
#' \itemize{
#'     \item{\strong{summary}: Summary of cross-validation, including mean and standard deviation of predictive ability (i.e. correlation between phenotype and GEBV, estimated on the validation fold, then averaged over replicates (nTimes), Time taken by the computation and number of markers}
#'     \item{\strong{cv}: Vector of predictive abilities averaged over nFolds, for each of the nTimes replicates}
#'     \item{\strong{sd}: Standard deviation of the nTimes predictive abilities}
#'     \item{\strong{MSEP}: Square root of the mean-squared error of prediction, averaged over Ntimes}
#'     \item{\strong{SDMSEP}: Standard deviation of the Square root of the mean-squared error of prediction, averaged over Ntimes} 
#'     \item{\strong{bv_table}: Matrix of dimension n x 4. Columns are: 
#'     \itemize{
#'         \item{Real BV, i.e. pheno vector}
#'         \item{Predict BV: the nx1 vector of GEBVs}
#'         \item{gpreSD: Standart deviation of estimated GEBV}
#'         \item{CD: coefficient of determination for each GEBV, estimated as sqrt Note that gpredSD and CD are only available for methods using the BGLR library, namely GBLUP, EGBLUP, BA,BB,BC,BL,RKHS and MKRKHS. These two columns contain NA for methods RF, RR, LASSO, EN and SVM.}
#'    }
#'  }
#' 
#' }
#' @examples
#' \donttest{
#' data(inra)
#' # Cross validation using GBLUP method
#' cv_gblup <- bwgs.cv(TRAIN47K, YieldBLUE, 
#'                        geno.impute.method = "mni",
#'                        predict.method = "gblup", 
#'                        nFolds = 10, 
#'                        nTimes = 1)
#' }
#' @export
bwgs.cv <- function(geno,pheno, FIXED = "NULL", 
                    MAXNA = 0.2,
                    MAF = 0.05,
                    pop.reduct.method="NULL",
                    sample.pop.size="NULL",
                    geno.reduct.method="NULL",
                    reduct.marker.size="NULL",
                    pval="NULL",
                    r2="NULL",
                    MAP="NULL",
                    geno.impute.method="NULL",
                    predict.method="NULL",
                    nFolds, 
                    nTimes)
{
  #(c)2015 louis.gautier.tran@gmail.com & gilles.charmet@clermont.inra.fr
  message("2017 BWGS  - Version 1.10.0 Release date: 31/10/2017")
  #message("2015 Gilles Charmet & Louis Gautier Tran)
  
  
  #message("")
  start.time <- Sys.time()
  message("Start time:")
  print(start.time)
  message("")
  pop.reduct.method <- toupper(pop.reduct.method)
  geno.reduct.method <- toupper(geno.reduct.method)
  
  
  reduct.size <- as.numeric(reduct.marker.size)
  geno.impute.method <- toupper(geno.impute.method)
  predict.method <- toupper(predict.method)
  #r2 is for LD reduct method
  #r2 and n.parts are for GMSLD or GMSLDV reduct method
  
  #if geno.impute.method ="ALL" and predict.method="ALL": bwgs.cv() will compare all methods
  
  #///////////////////////////////////////////////////////////////
  #STEP 0: select common lines in geno and pheno matrices
  # FILTER according to MAF and MAXNA
  #//////////////////////////////////////////////////////////////
  
  if(MAP!="NULL")  {
    MAPPED_markers=intersect(rownames(MAP),colnames(geno))
    if (length(MAPPED_markers) == 0) {
      stop("The row names of MAP is not matched with columns of geno")
    }
    MAP=MAP[MAPPED_markers,]
    geno=geno[,MAPPED_markers]
  }
  
  if(FIXED!="NULL") {
    genoFIX=intersect(rownames(FIXED),rownames(geno))
    
    FIXED=FIXED[genoFIX,]
    geno=geno[genoFIX,]
  }
  
  
  
  listGENO=rownames(geno)
  listPHENO=names(pheno)
  LIST=intersect(listGENO,listPHENO)
  
  if (length(LIST)==0)
  {stop("NO COMMON LINE BETWENN geno AND pheno")}
  if(length(LIST)>0)
  {
    geno=geno[LIST,]
    pheno=pheno[LIST]
    if(FIXED!="NULL"){FIXED=FIXED[LIST,]}
    
    new.pheno.size=length(LIST)
    message("Number of common lines between geno and pheno")
    print(new.pheno.size)
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
  
  if (pop.reduct.method=="NULL"|(toupper(as.character(pop.reduct.method))=="NULL")) {   
    rps <-0  
    message("random pop size not applied")
    print(rps)
    message("")
    
  } 
  
  
  
  if (pop.reduct.method=="RANDOM"|(toupper(as.character(pop.reduct.method))=="RANDOM")) {
    rps <- sample.pop.size
    
    message("random pop size=")
    print(rps)
    message("")
  }
  
  if (pop.reduct.method=="OPTI"|(toupper(as.character(pop.reduct.method))=="OPTI"))
  {
    #shrink_size <- as.numeric(reduct.size)
    rps <- as.numeric(sample.pop.size)
    
    
    message("optimized pop size=")
    print(rps)
    message("")	
  }
  
  
  message("POPULATION SAMPLING METHOD")
  print(pop.reduct.method)
  message("")
  
  
  
  
  
  #//////////////////////////////////////////////////////////////
  #Step 1: Imputation "MNI" or "EMI"
  #//////////////////////////////////////////////////////////////
  
  if((geno.impute.method=="NULL")|(geno.impute.method=="MNI")|(geno.impute.method=="EMI")) {
    
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
    if (FIXED!="NULL") {FIXED=round(MNI(FIXED))}
    
    
  } else {
    stop("Please choose an impute method: NULL, MNI, EMI ")
  } #If not chosen an impute method
  
  
  #///////////////////////////////////////////////////////////////
  #STEP 2: reduction OF THE DATA (reduire les donnees)
  #//////////////////////////////////////////////////////////////
  
  if((geno.reduct.method=="NULL")|(geno.reduct.method=="RMR")|(geno.reduct.method=="ANO")|(geno.reduct.method=="LD")|(geno.reduct.method=="ANO+LD")) { 
    
    if(geno.reduct.method=="NULL"){ 
      geno_shrink <- geno_impute 
      message("No reduction for genomic data.")
      new.geno.size <- dim(geno_shrink)
      
    }
    
    
    if(geno.reduct.method=="RMR"){   
      
      #if(is.null(reduct.size)) {stop("Please choose the size of columns for new genotypic data.")}
      
      if(!is.null(reduct.size)){
        #shrink_size <- as.numeric(reduct.size)
        time.rmr = system.time(geno_shrink <- RMR(geno_impute,reduct.size))
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
        
        time.ano = system.time(geno_shrink <- ANO(pheno,geno_impute,pval))
        time.ano.reduct = as.numeric(round(time.ano[3]/60,digits=2))
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
        
        time.ano = system.time(geno_shrinkANO <- ANO(pheno,geno_impute,pval))
        time.ano.reduct = as.numeric(round(time.ano[3]/60,digits=2))
        geno_shrinkANO=geno_shrinkANO[,!is.na(colnames(geno_shrinkANO))]
        
      }  else {stop("Please choose the p value for new genotypic data.")}
      
      if(MAP=="NULL")
      {
        stop("Please choose the r2 and/or MAP for LD reduction.")
      } # END OF LD reduction
      
      
      
      if(MAP!="NULL")
      {
        if (!is.null(r2))
        {
          MAP2=MAP[colnames(geno_shrinkANO),]
          time.chrld = system.time(geno_shrink <- CHROMLD(geno_shrinkANO,R2seuil=r2,MAP=MAP2))
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
        
        else {stop("Please choose the r2 and/or MAP for LD reduction.")}
      } # END OF LD reduction
      
    }
    
    if(geno.reduct.method=="LD")
      
    {
      if(MAP=="NULL")
      {
        
        {stop("Please choose a MAP for LD reduction.")}
      } # END OF LD reduction
      
      
      MAP2=MAP[colnames(geno),]
      
      time.chrld = system.time(geno_shrink <- CHROMLD(geno_impute,R2seuil=r2,MAP=MAP2))
      time.chrld.reduct = as.numeric(round(time.chrld[3]/60,digits=2))
      new.geno.size <- dim(geno_shrink)
      #geno_shrink
      message("Reduced by CHROMLD. New genotypic data dimension:")
      print(new.geno.size)
      message("")
      message("Time of reduction by CHROMLD (mins):")
      print(time.chrld.reduct)
      
      
      
      
    } # END OF LD reduction
    
    
    
    
  } # End if not chosen a reduct method
  
  else {stop("BWGS Warning! Please choose a valid reduct method. Ex.: RMR, AM, LD, NULL.")} #If not choose a reduct method
  
  #End of dimension reduction for the genomic matrix 
  
  
  
  
  #/////////////////////////////////////////////////////
  #STEP 3: CROSS VALIDATION - BEST MODEL IDENTIFICATION
  #/////////////////////////////////////////////////////
  if((predict.method=="EN")|(predict.method=="SVM")|(predict.method=="EGBLUP")|(predict.method=="RR")|(predict.method=="BA")|
     (predict.method=="BB")|(predict.method=="LASSO")|(predict.method=="BL")|(predict.method=="BC")|(predict.method=="GBLUP")|(predict.method=="BRNN")|
     (predict.method=="MKRKHS")|(predict.method=="RF")|(predict.method=="BRR")|(predict.method=="RKHS")|(predict.method=="ALL"))
  {
    
    
    
    if(predict.method=="GBLUP"){ 
      
      message("Predict by GBLUP...")
      # time.gblup = system.time(GBLUP <-runCV(pheno,geno_shrink,pop.reduct.method,rps,predictor_GBLUP,nFolds,nTimes))
      time.gblup = system.time(GBLUP <-runCV(pheno,geno_shrink,FIXED,pop.reduct.method,rps,predict_GBLUP,nFolds,nTimes))
      message("Predict by GBLUP...ended.")
      time.cal <- time.gblup
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(GBLUP,time.cal,ncol_geno_shrink)
      
    } 
    
    if(predict.method=="EGBLUP"){ 
      
      message("Predict by EGBLUP...")
      time.egblup = system.time(EGBLUP <-runCV(pheno,geno_shrink,FIXED,pop.reduct.method,rps=0,predict_EGBLUP,nFolds,nTimes))
      message("Predict by EGBLUP...ended.")
      time.cal <- time.egblup
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(EGBLUP,time.cal,ncol_geno_shrink)  
    }
    
    if(predict.method=="RF"){ 
      message("Predicting by RF...")
      message("Warning! The processing time for Random Forest Model is very long...")
      time.rf = system.time(RF <-runCV(pheno,geno_shrink,FIXED="NULL",pop.reduct.method,rps,predict_RF,nFolds,nTimes))
      message("Predict by RF...ended.")
      time.cal <- time.rf
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(RF,time.cal,ncol_geno_shrink)
      
    } 
    
    if(predict.method=="BRNN"){ 
      
      message("Predict by BRNN...")
      time.brnn = system.time(BRNN <-runCV(pheno,geno_shrink,FIXED="NULL",pop.reduct.method,rps,predict_BRNN,nFolds,nTimes))
      message("Predict by BRNN...ended.")
      time.cal <- time.brnn
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(BRNN,time.cal,ncol_geno_shrink)
      
    }
    
    if(predict.method=="RKHS"){ 
      message("Predicting by RKHS...")
      time.rkhs = system.time(RKHS <-runCV(pheno,geno_shrink,FIXED ="NULL",pop.reduct.method,rps,predict_RKHS,nFolds,nTimes))
      message("Predict by RKHS...ended.")
      time.cal <- time.rkhs
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(RKHS,time.cal,ncol_geno_shrink)
      
    } 
    
    if(predict.method=="MKRKHS"){ 
      message("Predicting by MKRKHS...")
      time.mkrkhs = system.time(MKRKHS <-runCV(pheno,geno_shrink,FIXED="NULL",pop.reduct.method,rps,predict_MKRKHS,nFolds,nTimes))
      message("Predict by MKRKHS...ended.")
      time.cal <- time.mkrkhs
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(MKRKHS,time.cal,ncol_geno_shrink)
      
    } 
    
    
    if(predict.method=="RR"){  #Ridge Regression#
      message("Predicting by RR (Ridge Regression)...")
      time.rr = system.time(RR <-runCV(pheno,geno_shrink,FIXED="NULL",pop.reduct.method,rps,predict_RR,nFolds,nTimes))
      message("Predict by RR...ended.")
      time.cal <- time.rr
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(RR,time.cal,ncol_geno_shrink)
      
    } 
    
    if(predict.method=="BRR"){  #Bayesian Ridge Regression#
      message("Predicting by BRR (Bayesian Ridge Regression)...")
      time.brr = system.time(BRR <-runCV(pheno,geno_shrink,FIXED,pop.reduct.method,rps,predict_BRR,nFolds,nTimes))
      message("Predict by BRR...ended.")
      time.cal <- time.brr
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(BRR,time.cal,ncol_geno_shrink)
      
    } 
    
    
    if(predict.method=="LASSO"){  #LASSO Regression#
      message("Predicting by LASSO ...")
      time.lasso= system.time(LASSO <-runCV(pheno,geno_shrink,FIXED="NULL",pop.reduct.method,rps,predict_Lasso,nFolds,nTimes))
      message("Predict by LASSO...ended.")
      time.cal <- time.lasso
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(LASSO,time.cal,ncol_geno_shrink)
    } 
    
    if(predict.method=="BL"){  #BLR regression#
      message("Predicting by Bayesian Lasso ...")
      time.bl= system.time(BL <-runCV(pheno,geno_shrink,FIXED,pop.reduct.method,rps,predict_BL,nFolds,nTimes))
      message("Predict by Bayesian Lasso...ended.")
      time.cal <- time.bl
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(BL,time.cal,ncol_geno_shrink)
    } 
    
    if(predict.method=="BA"){  #BayesianC regression#
      message("Predicting by BayesA ...")
      # time.ba= system.time(BA <-runCV(pheno,geno_shrink,pop.reduct.method,rps,predictor_BA,nFolds,nTimes))
      time.ba= system.time(BA <-runCV(pheno,geno_shrink,FIXED,pop.reduct.method,rps,predict_BA,nFolds,nTimes))
      message("Predict by BayesA...ended.")
      time.cal <- time.ba
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(BA,time.cal,ncol_geno_shrink)
    } 
    
    if(predict.method=="BB"){  #BayesB regression#
      message("Predicting by BayesB ...")
      time.bb= system.time(BB <-runCV(pheno,geno_shrink,FIXED,pop.reduct.method,rps,predict_BB,nFolds,nTimes))
      message("Predict by BayesB...ended.")
      time.cal <- time.bb
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(BB,time.cal,ncol_geno_shrink)
    } 
    
    if(predict.method=="BC"){  #BayesianC regression#
      message("Predicting by BayesC ...")
      time.bc= system.time(BC <-runCV(pheno,geno_shrink,FIXED,pop.reduct.method,rps,predict_BC,nFolds,nTimes))
      message("Predict by BayesC...ended.")
      time.cal <- time.bc
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(BC,time.cal,ncol_geno_shrink)
    } 
    
    if(predict.method=="EN"){  #Elastic-Net#
      message("Predicting by Elastic-Net ...")
      time.en = system.time(ElasticNet <-runCV(pheno,geno_shrink,FIXED="NULL",pop.reduct.method,rps,predict_ElasticNet,nFolds,nTimes))
      message("Predict by Elastic-Net...ended.")
      time.cal <- time.en
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(ElasticNet,time.cal,ncol_geno_shrink)
      
      
    } 
    
    if(predict.method=="SVM"){  #Support Vector Machine#
      message("Predicting by SVM ...")
      time.svm = system.time(SVM <-runCV(pheno,geno_shrink,FIXED="NULL",pop.reduct.method,rps,predict_SVM,nFolds,nTimes))
      message("Predict by SVM...ended.")
      time.cal <- time.svm
      ncol_geno_shrink <- ncol(geno_shrink)
      Results <- transfer(SVM,time.cal,ncol_geno_shrink)
    } 
    
    #FOR ALL PREDICTION METHODS:
    
    
    if (predict.method=="ALL")
    {
      
      message("Predict by all methods...")
      time.all = system.time(GSALL <-Run_All_Cross_Validation(pheno,geno_shrink,nFolds,nTimes))
      message("Predict by ALL METHODS...ended.")
      time.cal <- time.all
      ncol_geno_shrink <- ncol(geno_shrink)
      # Results <- transfer(GSALL,time.cal,ncol_geno_shrink)
      Results <- GSALL
      Results
    }
  } 
  
  #  else {stop("Please choose a predict method: EN,SVM,RR,EGBLUP, BA, BB, BC, BL, LASSO,GBLUP, RKHS, RF.")} # End of all predict methods
  
  else {stop("Please choose a predict method: EN,SVM,RR,LASSO,rrBLUP, RKHS, RF.")} # End of all predict methods
  
  
  Results
  
} 

