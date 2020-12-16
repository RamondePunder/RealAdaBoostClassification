#' Real Adaboost
#' 
#' Real Adaboost Classification based on Friedman, Hastie and Tibshirani (2000)
#' 
#' @author Ramon de Punder (University of Amsterdam, Tinbergen Institute)
#' @param mDataTrain matrix containig training data, on which weak learners are fitted.
#' @param mDataTest matrix containing test data, used for prediction and evaluation.
#' @param iM integer, total number of included weak learners, default is 10.
#' @param iDepth integer, maximal depth weak learners, default is 2.
#' @param iSeed integer, random seed used to do weighted sampling in each iteration step.
#' @param bVerbose boolean, prints weak learners when TRUE
#' @return A list containing class predictions, misclassification rate and a list including the base learners.
#' @export 
RealAdaBoostClass <- function(mDataTrain = NA, mDataTest = NA, iM = 10, iDepth = 2, iSeed = 1, bVerbose=TRUE){
  
  # Local functions [easier for package creation]
  Gini <- function(vY){
    "
    Purpose:
      Calculate Gini index for binary case
      
    Inputs:
      vY    vector, if vY is empty then return dGini=0
            Note: Consider binary predictor on which split has already been applied,
            then ChildLeft and ChildRight are vectors with only one value and therefore
            another split on those children based on the same regressor is pointless, since 
            then one of the grandchildren [Right one] of e.g. Childleft is empty. 
      
    Return:
      Gini index
    "
    dP1 <- sum(vY)/length(vY)
    dP0 <- 1- dP1
    dGini <- 1-dP1^2-dP0^2
    dGini <- replace(dGini, is.na(dGini), 0)
    return(dGini)
  }
  
  GiniImprov <- function(mParent=NA, sRegrName=NA, dSplit=NA){
    "
    Purpose:
      Calculate weighted improvement Gini when splitting sRegrName at dSplit
      
    Inputs:
      mParent       matrix, data before split
      sRegrName     string, regressor name
      dSplit        split value sRegrName
      
    Return:
      dDelta Gini   weighted improvement Gini index
    "
    
    # Initialisation
    vY <- mParent[,1]
    vX <- mParent[,sRegrName]
    bChildLeft <- vX < dSplit
    
    # Gini index parent and children
    dGiniParent <- Gini(vY)
    dGiniChildLeft <- Gini(vY[bChildLeft])
    dGiniChildRight <- Gini(vY[!bChildLeft])
    
    # Calculater weighted improvement
    dDeltaGini <- dGiniParent - sum(bChildLeft)/length(vY) * dGiniChildLeft  -  sum(!bChildLeft)/length(vY) * dGiniChildRight
    
    return(dDeltaGini)
  }
  
  OptimalSplitGini <- function(mParent=NA, lSplits=NA){
    "
    Purpose:
      Calculate optimal split in specific node of the tree
      
    Inputs:
      mParent     matrix, data parent node with dependent variable on first column
      lSplits     list, split values regressors
      
    Return: 
      lSplitOpt   list, optimal split: regressor name, treshold and Gini improvement
    "
    
    lRegrNames <- colnames(mParent)[-1]
    mSplit <- matrix(data=NA, nrow=length(lRegrNames), ncol=2)
    for(i in 1:length(lRegrNames)){
      sRegrName <- lRegrNames[i]
      vSplits <- lSplits[[sRegrName]]
      vGiniImpr <- matrix(data=NA, nrow=length(vSplits))
      for(j in 1:length(vSplits)){
        vGiniImpr[j,] <- GiniImprov(mParent, sRegrName, vSplits[j])
      }
      mSplit[i,] <- c(vSplits[which.max(vGiniImpr)], vGiniImpr[which.max(vGiniImpr),])
    }
    iOpt <- which.max(mSplit[,2])
    sRegrOpt <- lRegrNames[iOpt]
    iCat <- length(lSplits[[sRegrOpt]])+1
    iCount <- nrow(mParent)
    
    # Return: split variable name, count, ncat, index, improvement
    return(c(lRegrNames[iOpt], iCount, iCat, round(mSplit[iOpt,], 4)))
  }
  
  ClassificationTree <- function(mData=NA, iMinObs=10, iMaxDepth=2){
    "
    Purpose:
      Build tree based on relative Gini improvements
      
    Inputs:
      mData         matrix, data containing:
                      vY  vector, dependent variable on first column
                      mX  matrix, predictor variables
      iMinObs       integer, minimal number of observations before split
      iMaxDepth     integer, maximal depth tree
      
    Return:
      Split Table
    "
    # Initialisation
    mX <- mData[,-1]
    sXNames <- colnames(mX)
    iP <- length(sXNames)
    vNodes <- 1:(2^iMaxDepth - 1)
    iNNodes <- length(vNodes)
    vObsNodes <- rep(0, size=iNNodes) # number of obs. node i
    k <- 0 # level counter
    lNames <- list('node', 'level', 'variable', 'count,', 'ncat', 'index', 'improvement', 'p1_hat', 'y_hat')
    mSplitTable <- matrix(data=NA, nrow=iNNodes, ncol=length(lNames), dimnames=list(1:iNNodes,lNames)) 
    mSplitTable[,1] <- vNodes
    
    lSplits <- new.env() # like Python dictionary
    lDataNode <- list()  # data set each node
    lDataNode[[1]] <- mData 
    for(s in 1:iP){
      sX <- sXNames[s]
      vX_unique <- sort(unique(mX[,sX]))
      lSplits[[sX]] <- (vX_unique[-length(vX_unique)] + vX_unique[-1])/2
      if(length(vX_unique)==1) lSplits[[sX]] <- vX_unique # overwrite when vX only has one unique value
    }
    while(k < iMaxDepth){
      k <- k + 1
      for(i in 2^(k-1):(2^k -1)){ # loop through nodes at level k
        if(i>1){if(mSplitTable[floor(i/2),'variable'] == 'leaf'){next}
        }
        if(max(nrow(lDataNode[[i]]), 0) >= iMinObs){ 
          #note: max for cases where iMinObs is large and nrow(lDataNode[[i]] may become NULL
          
          mParent= lDataNode[[i]]
          
          # Update Split Table
          mSplitTable[i,] <- c(i, k, OptimalSplitGini(mParent= mParent, lSplits=lSplits),NA,NA)
          
          # Construct Children
          bChildLeft <- mParent[,mSplitTable[i,'variable']] < as.numeric(mSplitTable[i,'index'])
          lDataNode[[2 * i]] <- mParent[bChildLeft,]
          lDataNode[[2 * i + 1]] <- mParent[!bChildLeft,]
        }
        else{
          mSplitTable[i,] <- c(i, k, 'leaf', max(nrow(lDataNode[[i]]), 0), NA, NA, NA,NA,NA)
        }
      } 
    }
    # Identify nodes which have leaf as parent
    bLeafParent <- !is.na(mSplitTable[,'level'])
    
    # All nodes at at highest level which are not yet set to leaves by iMinObs criterion are leaves as well:
    vLeaves <- 2^(iMaxDepth-1):(2^iMaxDepth -1) 
    for(i in 1:length(vLeaves)){
      n <- vLeaves[i]
      iCount <- if(bLeafParent[n]==TRUE){
        iCount <- as.integer(length(lDataNode[[n]])!=0) * max(nrow(lDataNode[[n]]), 1)  
        # max for leaves with single value and return 0 for empty leaves
      }
      else{
        iCount <- NA
      }
      mSplitTable[n,] <- c(n, iMaxDepth, 'leaf', iCount, NA, NA, NA,NA,NA)
    } 
    
    mSplitTableRed <- mSplitTable[bLeafParent,] # drop nodes with leaf as parent
    for(i in 1:length(mSplitTableRed[,'node'])){
      n <- as.integer(mSplitTableRed[,'node'][i])
      mDataNode <- matrix(lDataNode[[n]], ncol=ncol(mData)) # matrix to make robust against vector object lDataNode[[n]]
      mSplitTableRed[i,ncol(mSplitTableRed)-1] <-  round(sum(pmax(mDataNode[,1], 0))/length(mDataNode[,1]),4)  # proportion of ones
      dIndUpper <- round(as.numeric(mSplitTableRed[i,ncol(mSplitTableRed)-1],0))
      mSplitTableRed[i,ncol(mSplitTableRed)] <- dIndUpper * 1 + (1-dIndUpper) * -1
      if(length(mDataNode[,1]) == 0) mSplitTableRed[i,c(ncol(mSplitTableRed)-1, ncol(mSplitTableRed))] = c(1/2, NA)
    }
    return(mSplitTableRed)
  }
  
  ClassificationTreePred <- function(mSplitTable=NA, vData=NA, bProb= TRUE){
    "
    Purpose:
      Prediction and error calculation based on single row data set
      
    Inputs:
      mSplitTable   matrix, split table fitted tree
      vData         vector, single row data set
      bProb         boolean, probability returns class probability when TRUE
      
    Return:
      dYhat         prediction
      dError        prediction error
    "
    n <- 1
    while(n <= max(mSplitTable[,'node'])){
      vSplitTableSelect <- mSplitTable[mSplitTable[,'node']==n,]
      if(vSplitTableSelect['variable']=='leaf'){
        dYhat= as.numeric(vSplitTableSelect['y_hat']) # prediction is majority vote in leaf
        dPhat= as.numeric(vSplitTableSelect['p1_hat'])
        break
      }
      else{
        if(vData[vSplitTableSelect['variable']] < vSplitTableSelect['index']){
          n <- 2*n    # Left child
        }
        else{
          n <- 2*n +1 # Right child
        }
      }
    }
    dError <- vData[1]-dYhat
    return(list('y_hat'=dYhat, 'p_1hat'=dPhat, 'error'=dError))}
  
  
  
  
  iN <- nrow(mDataTrain)
  set.seed(iSeed)
  vWeights <- rep(1/iN,iN)
  vPred <- matrix(0, nrow = nrow(mDataTest))
  vF <- matrix(NA, nrow = iM)
  lWeakLearners <- list()
  for(m in 1:iM){
    
    # Obtain weighted sample from training sample
    vIdxWeighted <- sample(iN, iN, prob=vWeights, replace = TRUE)
    mDataTrain_m <- mDataTrain[vIdxWeighted, ]
    
    # Classifier based on training data and current weights
    mSplitTable <- ClassificationTree(mDataTrain_m, iMinObs=2, iMaxDepth=iDepth)
    lWeakLearners[[m]] <- mSplitTable
    
    Fm <- function(vX){
      # Set precision, to avoid issues in log(p/(1-p)) for p=1 and p=0 [see original paper Friedman et al (2000)]
      dEps = 0.001
      
      # Note: we can let the corresponding y be 1, as we only use the p_1hat [and not the error]
      # called by ClassificationTreePred
      vData <- setNames(c(1, vX), colnames(mDataTrain))
      dPm <- min(max(ClassificationTreePred(mSplitTable=mSplitTable, vData=vData)$p_1hat, dEps), 1-dEps)
      dFm <- 1/2  * log(dPm / (1 - dPm))
      return(dFm)
    }
    
    # Update weights
    vWeights <- vWeights * exp(-mDataTrain_m[,1] * apply(mDataTrain_m[,-1], 1, Fm))
    
    # Calculate output classifier part
    vPred <- vPred + apply(mDataTest[,-1], 1, Fm)
  }
  if(bVerbose == TRUE) print(lWeakLearners)
  vPred <- sign(vPred)
  dError <- 1 - sum(vPred==mDataTest[, 1])/length(vPred) # misclassification rate
  return(list('Pred'=vPred, 'Error'=dError, 'Learners'=lWeakLearners))
}

#' Real Adaboost Cross Validation
#' 
#' Cross Validation for Real Adaboost Classification based on Friedman, Hastie and Tibshirani (2000)
#' @param lCVData list, has length equal to number of folds, each list element lists a train and test matrix of data 
#' @param lArgs list, contains other parameters needed for RealAdaBoostClass than to be tuned by CV: Depth and Seed.
#' @param vM vector, grid for M, the total number of included weak learners
#' @param sName, string, characteristic part file name
#' @param bPlot, boolean, creates plot and saves tikzfile in folder Figures when true.
#' @return CV value for M, optionally together with corresponding standard error when bPlot is FALSE.
#' @export
RealAdaBoostClass_CV <- function(lCVData = NA, lArgs = NA, vM= seq(1,15), sName ='', bPlot = TRUE){

  # Packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(rstudioapi,    # set working directory
                 doParallel,    # detect cores
                 ggplot2,       # plotting
                 reshape2,       # melt function
                 tikzDevice)    # create tikzfiles 
  
  # Initialisation
  iFold <- length(lCVData)
  iCores <- detectCores()
  iN <- nrow(lCVData[[1]][[1]]) + nrow(lCVData[[1]][[2]])
  iDepth <- lArgs[['Depth']]
  iSeed <- lArgs[['Seed']]
  
  # Define function handle
  Error <- function(iM){
    vError <- matrix(NA, nrow = iFold)
    for(i in 1: iFold){
      # Use misclassification rate as performance measure
      vError[i] <- RealAdaBoostClass(mDataTrain = lCVData[[i]]$mDataTrain, mDataTest = lCVData[[i]]$mDataTest, iM = iM, iDepth = iDepth, iSeed = iSeed, bVerbose=FALSE)$Error
    }
    dError <- mean(vError) # overall misclassification rate is the mean of the misclassification rates of the folds
    dSE <- sqrt((dError * (1-dError))/ iN)   # variance is p(1-p)/n, use phat to operationalise
    return(c(dError, dSE))
  }
  
  # Run RealAdaBoost in parallel
  mError <- matrix(unlist(mclapply(vM, Error, mc.cores = iCores)), nrow = length(vM), ncol = 2, byrow = TRUE)
  mCVFig <- cbind(vM, mError)
  
  # RMSE minimising lambda on grid
  iMMin <- mCVFig[mCVFig[, 2] == min(mCVFig[, 2])][1] 
  
  # Standard deviation rule
  dUpper <- mCVFig[mCVFig[, 1] == iMMin][2] + mCVFig[mCVFig[, 1] == iMMin][3]
  iMOpt <- min(mCVFig[mCVFig[, 2] < dUpper, 1]) # note: min for least complex model
  
  if(bPlot == TRUE){
    cat('Optimal value', sprintf(": %i", iFold), '-fold CV: ' ,'M',  sprintf(": %.4f.", iMOpt), "\n", sep = "")
    
    # Generate figure with error bars
    dfDataFig <- as.data.frame(mCVFig)
    colnames(dfDataFig) = c('M', 'Error', 'sd')
    setwd(paste(dirname(getActiveDocumentContext()$path), '/Figures', sep =''))
    tikz(file = paste("CV_M_", sName, sep = ''), width = 6,height = 5)
    pPlot <- ggplot(dfDataFig, aes_string(x = 'M', y = 'Error'), colour = variable, group = variable) +
      geom_errorbar(aes(ymin = Error - sd, ymax = Error + sd), width = .2, position = position_dodge(0.05), colour = 'grey') +
      geom_line(colour= 'blue') + geom_point(colour= 'blue') + theme_classic() + #scale_x_log10() +
      geom_vline(xintercept = iMMin, linetype="dashed", color = "black", size=0.5) +
      geom_vline(xintercept = iMOpt, linetype="dashed", color = "black", size=0.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size=16, color = 'black'),axis.title = element_text(size=16))
    print(pPlot)
    dev.off()
    print(pPlot)
    return(iMOpt)
  }
  else{
    return(c(iMOpt, mCVFig[mCVFig[,1]==iMOpt,2]))
  }
}


