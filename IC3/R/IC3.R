#' Title
#'
#' @param A  Gene expression count matrix, each row represents a cell and each colume represents a gene.
#' @param cellinfo  Cell information matrix, each row represents a cell.The first column represents the x-axis coordinates of the cell, and the second column represents the y-axis coordinates of the cell. The third column represents the cell type of the cell. The row name is the name of the cell, and its order should be consistent with the order of gene expression matrix A.
#' @param lrinfo lrinfo: Ligand-Receptor Database, each row represents a ligand-receptor pair. Ligand is in the first column.
#' @param alpha   Distance proportional threshold, cells with a distance exceeding the threshold are considered to have no interaction with each other. Default value is 0.02.
#' @param selfinter  Consider self communication or not, default is TRUE.
#' @param minitr The minimum iteration number. Default is 10.
#' @param maxitr The maximum iteration number. Default is 100.
#' @param minbeta The minimum value of parameter beta. Default is 0.
#' @return The first term: the communication probability with cell type level; The second term: the communication probability with single cell level. The 3-5th term: parameter estimation of lambda;beta;r=(r0,r1,r2). The 6th term: Q function value after each iteration
#' @export
#'
#' @examples IC3(A, cellinfo, lrinfo)
IC3 <- function(A, cellinfo, lrinfo, alpha = 0.02, selfinter = TRUE, minitr = 10, maxitr = 100, minbeta = 0) {
  library(progress)
  library(stringr)
  A <- as.matrix(A)
  if (length(which(A %% 1 > 0)) > 0) {
    return("error:there are non-intergers in your count matrix")
  }
  if (length(which(A < 0)) > 0) {
    return("error:there are negative numbers in your count matrix")
  }
  cellnum <- dim(A)[1]
  genenum <- dim(A)[2]
  if (dim(cellinfo)[1] != cellnum) {
    return("error:the first dimension of cell information matrix is not equal to the first dimension of A")
  }
  if (dim(cellinfo)[2] != 3) {
    return("error:the second dimension of cell information matrix is not equal to 3")
  }
  if (dim(lrinfo)[2] != 2) {
    return("error:the second dimension of ligand-receptor matrix is not equal to 2")
  }
  text <- paste("Data uploaded successfully, there are ",cellnum," cells and",genenum," genes."," The ligand receptor database has ",dim(lrinfo)[1]," different pairs.")
  print(text)
  location <- data.frame(cellinfo[,1:2])
  location[,1] <- as.numeric(location[,1]);
  location[,2] <- as.numeric(location[,2]);
  location <- as.matrix(location)
  cellname <- rownames(cellinfo)
  pairnum <- alpha*cellnum*cellnum/2
  candidatepair <- matrix(100000,pairnum,3);
  nearpair <- matrix(0,cellnum,3)
  pb <- progress_bar$new(total = cellnum)
  print("Calculating the distances between cells")
  for (i in 1:cellnum)
  {
      pb$tick()
      dian <- as.numeric(location[i,])
      dians <- as.matrix(location)
      distances <- ((dians[,1] - dian[1])^2+ (dians[,2] - dian[2])^2)^(0.5)
      nearpair[i,1] <- i;
      nearpair[i,2] <- order(distances)[2];
      nearpair[i,3] <- sort(distances)[2];
      thispair <- matrix(0,cellnum,3);
      thispair[,1] <- rep(i,cellnum);
      thispair[,2] <- 1:cellnum; 
      thispair[,3] <- distances;
      thispair <- thispair[-i,];
      candidatepair <- rbind(candidatepair,thispair)
      candidatepair <- candidatepair[order(candidatepair[,3])[1:pairnum],]
  }
  allpair <- rbind(candidatepair,nearpair)
  threshold <-  max(candidatepair[,3]);
  for (i in 1:nrow(allpair)) {
    if (allpair[i, 2] < allpair[i, 1]) {
      temp <- allpair[i, 1]
      allpair[i, 1] <- allpair[i, 2]
      allpair[i, 2] <- temp
    } 
  }  
  duplicate_rows <- duplicated(allpair)
  allpair_unique <- allpair[!duplicate_rows, ]
  if (selfinter)
  {
    selfpair <- matrix(0,cellnum,3);
    selfpair[,1] <- 1:cellnum;selfpair[,2] <- 1:cellnum;
    cellpair <- rbind(allpair_unique,selfpair); 
  }
  if (selfinter == FALSE)
  {
     cellpair <- allpair_unique; 
  }
  cellpairnum <- dim(cellpair)[1]
  text <- paste("The interaction ratio is ",alpha,".The interaction distance threshold is ", threshold)
  print(text) 
  text <- paste("There are ",cellpairnum," possible interact cell pairs.")
  print(text)
  cellneighbor <- list()
  for (i in 1:cellnum)
  {
     neighborone <- cellpair[which(cellpair[,1]==i),2]
     neighbortwo <- cellpair[which(cellpair[,2]==i),1]  
     cellneighbor[[i]] <- union(neighborone,neighbortwo)
  }
  library(Matrix)
  cellpairindex  <- Matrix(data = 0L, nrow=cellnum, ncol = cellnum, sparse = TRUE)
  for (i in 1:cellpairnum)
  {
     cellpairindex[cellpair[i, 1], cellpair[i, 2]] <- i
     cellpairindex[cellpair[i, 2], cellpair[i, 1]] <- i
  }  
  cellname <- as.character(rownames(A))
  genename <- as.character(colnames(A))
  lrinfo <- as.matrix(lrinfo)
  s <- 0
  lrpair <- matrix(0, dim(lrinfo)[1], 2)
  for (i in 1:dim(lrinfo)[1])
  {
    if (lrinfo[i, 1] %in% genename && lrinfo[i, 2] %in% genename) {
      s <- s + 1
      lrpair[s, 1] <- which(genename == lrinfo[i, 1])
      lrpair[s, 2] <- which(genename == lrinfo[i, 2])
    }
  }
  lrpair <- lrpair[1:s, ]
  genepairnum <- s
  text <- paste("There are ", genepairnum, "ligand-receptor pairs in the gene list")
  print(text)
  allgene <- sort(unique(as.numeric(lrpair)))
  ligand <- sort(unique(as.numeric(lrpair[, 1])))
  receptor <- sort(unique(as.numeric(lrpair[, 2])))
  lrname <- genename[allgene]
  y <- A[, allgene]
  genenum <- dim(y)[2]
  cellnum <- dim(y)[1]
  colnames(y) <- lrname
  geneneighbor <- list()
  for (i in 1:length(allgene))
  {
    if (allgene[i] %in% ligand) {
      lrindex <- which(lrpair[, 1] == allgene[i])
      receptorindex <- lrpair[lrindex, 2]
      receptorname <- genename[receptorindex]
      geneneighbor[[i]] <- which(lrname %in% receptorname)
    }
    if (allgene[i] %in% receptor) {
      lrindex <- which(lrpair[, 2] == allgene[i])
      ligandindex <- lrpair[lrindex, 1]
      ligandname <- genename[ligandindex]
      geneneighbor[[i]] <- which(lrname %in% ligandname)
    }
  }
  genepair <- matrix(0, genepairnum, 2)
  for (i in 1:genepairnum)
  {
    genepair[i, 1] <- which(lrname == genename[lrpair[i, 1]])
    genepair[i, 2] <- which(lrname == genename[lrpair[i, 2]])
  }
  genepairindex <- matrix(0, dim(y)[2], dim(y)[2])
  for (i in 1:dim(genepair)[1])
  {
    genepairindex[genepair[i, 1], genepair[i, 2]] <- i
    genepairindex[genepair[i, 2], genepair[i, 1]] <- i
  }
  geneneighborindex <- list()
  for (i in 1:genenum)
  {
    geneneighborindex[[i]] <- genepairindex[i, geneneighbor[[i]]]
  }
  celltype <- cellinfo[, 3]
  typename <- unique(celltype)
  typeindex <- rep(0, cellnum)
  for (i in 1:cellnum)
  {
    typeindex[i] <- which(typename == celltype[i])
  }
  typenum <- length(typename)
  ddata <- cellpair[,3]
  numc <- rowMeans(y)
  for (i in 1:length(numc)) {
    numc[i] <- max(0.001, numc[i])
  }
  conll <- function(ythis, cell0, gene0, ethis) {
    type0 <- typeindex[cell0]
    result <- ythis * log(lambda[type0, gene0] * numc[cell0]) - log(factorial(ythis)) ## No minus lambda[type0,gene0]*numc[cell0] because we need to plus it when  calculate the normalization.
    for (c in cellneighbor[[cell0]])
    {
      if (ethis[cellpairindex[cell0, c]] != 0) {
        for (g in geneneighbor[[gene0]])
        {
          result <- result + ethis[cellpairindex[cell0, c]] * beta[genepairindex[gene0, g]] * atan(ythis) * atan(y[c, g])
        }
      }
    }
    return(result)
  }
  lambda <- matrix(0, typenum, genenum)
  chulambda <- matrix(0, typenum, genenum)
  for (i in 1:typenum)
  {
    for (j in 1:genenum)
    {
      chulambda[i, j] <- sum(y[which(typeindex == i), j]) / sum(numc[which(typeindex == i)])
      lambda[i, j] <- max(chulambda[i, j], 0.001)
    }
  }
  beta <- rep(minbeta, genepairnum)
  CInum <- matrix(0, typenum, typenum)
  for (i in 1:cellpairnum)
  {
    t1 <- min(typeindex[cellpair[i, 1]], typeindex[cellpair[i, 2]])
    t2 <- max(typeindex[cellpair[i, 1]], typeindex[cellpair[i, 2]])
    CInum[t1, t2] <- CInum[t1, t2] + 1
  }
  ECI <- matrix(0, typenum, typenum)
  for (i in 1:typenum)
  {
    for (j in i:typenum)
    {
      if (CInum[i, j] > 0) {
        ECI[i, j] <- 0.5
      }
    }
  }
  logit <- function(x) {
    y <- log(x) - log(1 - x)
    return(y)
  }
  r0 <- 0.01 * mean(ddata)
  r1 <- 1
  r2 <- -0.01
  e <- rep(0, cellpairnum)
  for (i in 1:cellpairnum)
  {
    t1 <- min(typeindex[cellpair[i, 1]], typeindex[cellpair[i, 2]])
    t2 <- max(typeindex[cellpair[i, 1]], typeindex[cellpair[i, 2]])
    e[i] <- exp(r0 + r1 * logit(ECI[t1, t2]) + r2 * ddata[i]) / (1 + exp(r0 + r1 * logit(ECI[t1, t2]) + r2 * ddata[i]))
  }
  hatW <- rep(0, cellpairnum)
  for (i in 1:cellpairnum) {
    if (e[i] > 0.5) {
      hatW[i] <- 1
    }
  }
  Q <- function(e, y, r0, r1, r2, ECI, W, lambda, beta) {
    logit <- function(x) {
      y <- log(x) - log(1 - x)
      return(y)
    }
    logl <- 0
    for (i in 1:cellpairnum)
    {
      cell1 <- cellpair[i, 1]
      cell2 <- cellpair[i, 2]
      type1 <- min(typeindex[cell1], typeindex[cell2])
      type2 <- max(typeindex[cell1], typeindex[cell2])
      logl <- logl + e[i] * (r0 + r1 * logit(ECI[type1, type2]) + r2 * ddata[i]) - log(1 + exp(r0 + r1 * logit(ECI[type1, type2]) + r2 * ddata[i]))
    }
    qconll <- function(ythis, cell0, gene0) {
      result <- ythis * log(lambda[typeindex[cell0], gene0] * numc[cell0]) - log(factorial(ythis))
      for (c in cellneighbor[[cell0]])
      {
        for (g in geneneighbor[[gene0]])
        {
          result <- result + e[cellpairindex[cell0, c]] * beta[genepairindex[gene0, g]] * atan(ythis) * atan(y[c, g])
        }
      }
      return(result)
    }
    for (cell in 1:cellnum)
    {
      for (gene in 1:genenum)
      {
        m <- max(y[cell, gene], 10)
        B <- qconll(0:m, cell0 = cell, gene0 = gene)
        logl <- logl + B[y[cell, gene] + 1] - log(sum(exp(B)))
      }
    }
    return(logl)
  }
  library(stringr)
  lljilu <- matrix(0, maxitr, 5)
  for (itr in 1:maxitr)
  {
    print(paste("The", itr, "iterate", sep = " "))
    print("ICM step")
    hatW <- rep(0, cellpairnum)
    for (i in 1:cellpairnum) {
      if (e[i] > 0.5) {
        hatW[i] <- 1
      }
    }
    print("E step")
    pb <- progress_bar$new(total = cellpairnum)
    logit <- function(x) {
      y <- log(x) - log(1 - x)
      return(y)
    }
    newe <- rep(0, cellpairnum)
    for (i in 1:cellpairnum)
    {
      pb$tick()
      cell1 <- cellpair[i, 1]
      cell2 <- cellpair[i, 2]
      type1 <- min(typeindex[cell1], typeindex[cell2])
      type2 <- max(typeindex[cell1], typeindex[cell2])
      cellnum <- dim(y)[1]
      genenum <- dim(y)[2]
      logb <- 0
      change <- rep(0, genenum)
      for (g in 1:genenum)
      {
        old <- logb
        for (ng in geneneighbor[[g]]) {
          logb <- logb + beta[genepairindex[g, ng]] * atan(y[cell1, g]) * atan(y[cell2, ng])
        }
        for (ng in geneneighbor[[g]]) {
          logb <- logb + beta[genepairindex[g, ng]] * atan(y[cell2, g]) * atan(y[cell1, ng])
        }
        zone <- hatW
        zone[i] <- 1
        zzero <- hatW
        zzero[i] <- 0
        logb <- logb + log(sum(exp(conll(0:10, cell1, g, zzero)))) + log(sum(exp(conll(0:10, cell2, g, zzero))))
        logb <- logb - log(sum(exp(conll(0:10, cell1, g, zone)))) - log(sum(exp(conll(0:10, cell2, g, zone))))
        change[g] <- logb - old
      }
      b <- exp(logb)
      b <- b * exp(r0 + r1 * logit(ECI[type1, type2]) + r2 * ddata[i])
      newe[i] <- b / (1 + b)
    }
    e <- newe
    print("M step")
    print("Maximize beta")
    daoll <- function(ythis, cell0, gene2, ethis) {
      result <- 0
      for (c in cellneighbor[[cell0]])
      {
        result <- result + ethis[cellpairindex[cell0, c]] * atan(ythis) * atan(y[c, gene2])
      }
      return(result)
    }
    newbeta <- beta
    pb <- progress_bar$new(total = genepairnum)
    for (i in 1:genepairnum)
    {
      pb$tick()
      gene1 <- genepair[i, 1]
      gene2 <- genepair[i, 2]
      daoyi <- 0
      daoer <- 0
      for (c in 1:cellnum)
      {
        for (nc in cellneighbor[[c]]) {
          daoyi <- daoyi + 2 * e[cellpairindex[c, nc]] * atan(y[c, gene1]) * atan(y[nc, gene2])
        }
        Zlist <- conll(0:10, c, gene1, e)
        addlist <- daoll(0:10, c, gene2, e)
        B0gene12 <- sum(exp(Zlist))
        B1gene12 <- sum(exp(Zlist) * addlist)
        B2gene12 <- sum(exp(Zlist) * addlist^2)
        Zlist <- conll(0:10, c, gene2, e)
        addlist <- daoll(0:10, c, gene1, e)
        B0gene21 <- sum(exp(Zlist))
        B1gene21 <- sum(exp(Zlist) * addlist)
        B2gene21 <- sum(exp(Zlist) * addlist^2)
        daoyi <- daoyi - B1gene12 / B0gene12 - B1gene21 / B0gene21
        daoer <- daoer - (B2gene12 * B0gene12 - (B1gene12)^2) / (B0gene12^2)
        daoer <- daoer - (B2gene21 * B0gene21 - (B1gene21)^2) / (B0gene21^2)
      }
      kk <- daoyi / daoer
      if (is.na(kk)) {
        kk <- 0
      }
      if (kk > 5) {
        kk <- 5
      }
      if (kk < (-5)) {
        kk <- -5
      }
      newbeta[i] <- max(minbeta, beta[i] - 0.5 * kk)
    }
    candidatesulv <- c(1,0.5,0.2,0.1,0);
    lljilu[itr,1] <- Q(e, y, r0, r1, r2, ECI, hatW, lambda, beta)
    yuanll <- lljilu[itr,1]
    for (tt in 1:5)
    {
       candidatebeta <- candidatesulv[tt] * newbeta + (1-candidatesulv[tt]) * beta
       nowll <- Q(e, y, r0, r1, r2, ECI, hatW, lambda, candidatebeta)
       if (nowll > yuanll)
       {
           break;
       }
    }
    lljilu[itr,2] <- nowll
    newbeta <- candidatebeta
    chabeta <- max(abs(beta - newbeta))
    beta <- newbeta
    print("Maximize lambda")
    newlambda <- lambda
    dconll <- function(ythis, cell0, gene0, ethis, k) {
      type0 <- typeindex[cell0]
      result <- ythis * log(lambda[type0, gene0]) - log(factorial(ythis))
      for (c in cellneighbor[[cell0]])
      {
        if (ethis[cellpairindex[cell0, c]] != 0) {
          for (g in geneneighbor[[gene0]])
          {
            result <- result + ethis[cellpairindex[cell0, c]] * beta[genepairindex[gene0, g]] * atan(ythis + k) * atan(y[c, g])
          }
        }
      }
      return(result)
    }
    pb <- progress_bar$new(total = genenum * typenum)
    for (t in 1:typenum)
    {
      for (g in 1:genenum)
      {
        pb$tick()
        if (chulambda[t, g] != 0) {
          daoyi <- 0
          daoer <- 0
          typecell <- which(typeindex == t)
          for (c in typecell)
          {
            A0 <- sum(exp(dconll(0:10, c, g, e, 0)))
            A1 <- sum(exp(dconll(0:10, c, g, e, 1)))
            A2 <- sum(exp(dconll(0:10, c, g, e, 2)))
            daoyi <- daoyi + y[c, g] / lambda[t, g] - numc[c] * A1 / A0
            daoer <- daoer - y[c, g] / lambda[t, g]^2 - numc[c]^2 * (A2 * A0 - A1^2) / A0^2
          }
          daoer <- -abs(daoer)
          kk <- daoyi / daoer
          if (is.na(kk)) {
            kk <- 0
          }
          if (kk > 5) {
            kk <- 5
          }
          if (kk < (-5)) {
            kk <- -5
          }
          newlambda[t, g] <- max(0.001, min(chulambda[t, g], (lambda[t, g] - 0.5 * kk)))
        }
      }
    }
    candidatesulv <- c(1,0.5,0.2,0.1,0);
    yuanll <- lljilu[itr,2];
    for (tt in 1:5)
    {
       candidatelambda <- candidatesulv[tt] * newlambda + (1-candidatesulv[tt]) * lambda
       nowll <- Q(e, y, r0, r1, r2, ECI, hatW, candidatelambda, beta)
       if (nowll > yuanll)
       {
           break;
       }
    }
    lljilu[itr,3] <- nowll;
    newlambda <- candidatelambda
    chalambda <- max(abs(newlambda - lambda))
    lambda <- newlambda
    print("Maximize I")
    daoyi <- matrix(0, typenum, typenum)
    daoer <- matrix(0, typenum, typenum)
    logit <- function(x) {
      y <- log(x) - log(1 - x)
      return(y)
    }
    pb <- progress_bar$new(total = cellpairnum)
    for (i in 1:cellpairnum)
    {
      pb$tick()
      cell1 <- cellpair[i, 1]
      cell2 <- cellpair[i, 2]
      type1 <- min(typeindex[cell1], typeindex[cell2])
      type2 <- max(typeindex[cell1], typeindex[cell2])
      jiben <- r0 + r1 * logit(ECI[type1, type2]) + r2 * ddata[i]
      daoyi[type1, type2] <- daoyi[type1, type2] + e[i] * r1 * (1 / ECI[type1, type2] + 1 / (1 - ECI[type1, type2]))
      fenzi <- exp(jiben) * r1 * (1 / ECI[type1, type2] + 1 / (1 - ECI[type1, type2]))
      fenmu <- 1 + exp(jiben)
      daoyi[type1, type2] <- daoyi[type1, type2] - fenzi / fenmu
      daoer[type1, type2] <- daoer[type1, type2] + r1 * e[i] * (-1 / (ECI[type1, type2])^2 + 1 / (1 - ECI[type1, type2])^2)
      fenzi <- exp(jiben) * r1^2 * (1 / ECI[type1, type2] + 1 / (1 - ECI[type1, type2]))^2
      fenmu <- (1 + exp(jiben))^2
      daoer[type1, type2] <- daoer[type1, type2] - fenzi / fenmu
      fenzi <- exp(jiben) * r1 * (-1 / (ECI[type1, type2])^2 + 1 / (1 - ECI[type1, type2])^2)
      fenmu <- 1 + exp(jiben)
      daoer[type1, type2] <- daoer[type1, type2] - fenzi / fenmu
    }
    newECI <- ECI
    for (i in 1:typenum)
    {
      for (j in 1:typenum)
      {
        if (CInum[i, j] > 0) {
          kk <- daoyi[i, j] / daoer[i, j]
          if (is.na(kk)) {
            kk <- 0
          }
          if (kk > 5) {
            kk <- 5
          }
          if (kk < (-5)) {
            kk <- -5
          }
          newECI[i, j] <- max(0.001, min(0.999, ECI[i, j] - 0.5 * kk))
        }
      }
    }
    candidatesulv <- c(1,0.5,0.2,0.1,0);
    yuanll <- lljilu[itr,3]
    for (tt in 1:5)
    {
       candidateECI <- ECI
       for (di in 1:typenum)
       {
          for (dj in 1:typenum)
          {
              candidateECI[di,dj]=candidatesulv[tt] * newECI[di,dj] + (1-candidatesulv[tt]) * ECI[di,dj]
          }
       }
       nowll <- Q(e, y, r0, r1, r2, candidateECI, hatW, lambda, beta)
       if (nowll > yuanll)
       {
           break;
       }
    }
    lljilu[itr,4] <- nowll;
    newECI <- candidateECI
    chaECI <- max(abs(newECI - ECI))
    ECI <- newECI
    print("Maximize r0,r1,r2")
    logit <- function(x) {
      y <- log(x) - log(1 - x)
      return(y)
    }
    pb <- progress_bar$new(total = cellpairnum)
    zerodaoyi <- 0
    zerodaoer <- 0
    onedaoyi <- 0
    onedaoer <- 0
    twodaoyi <- 0
    twodaoer <- 0
    for (i in 1:cellpairnum)
    {
      pb$tick()
      cell1 <- cellpair[i, 1]
      cell2 <- cellpair[i, 2]
      type1 <- min(typeindex[cell1], typeindex[cell2])
      type2 <- max(typeindex[cell1], typeindex[cell2])
      jiben <- r0 + r1 * logit(ECI[type1, type2]) + r2 * ddata[i]
      zerodaoyi <- zerodaoyi + e[i] - exp(jiben) / (1 + exp(jiben))
      zerodaoer <- zerodaoer - exp(jiben) / (1 + exp(jiben))^2
      onedaoyi <- onedaoyi + (e[i] - exp(jiben) / (1 + exp(jiben))) * logit(ECI[type1, type2])
      onedaoer <- onedaoer - logit(ECI[type1, type2])^2 * exp(jiben) / (1 + exp(jiben))^2
      twodaoyi <- twodaoyi + (e[i] - exp(jiben) / (1 + exp(jiben))) * ddata[i]
      twodaoer <- twodaoer - ddata[i]^2 * exp(jiben) / (1 + exp(jiben))^2
    }
    kk <- zerodaoyi / zerodaoer
    if (is.na(kk)) {
      kk <- 0
    }
    if (kk > 5) {
      kk <- 5
    }
    if (kk < (-5)) {
      kk <- -5
    }
    newr0 <- r0 - 0.5 * kk
    kk <- onedaoyi / onedaoer
    if (is.na(kk)) {
      kk <- 0
    }
    if (kk > 5) {
      kk <- 5
    }
    if (kk < (-5)) {
      kk <- -5
    }
    newr1 <- r1 - 0.5 * kk
    kk <- twodaoyi / twodaoer
    if (is.na(kk)) {
      kk <- 0
    }
    if (kk > 5) {
      kk <- 5
    }
    if (kk < (-5)) {
      kk <- -5
    }
    newr2 <- r2 - 0.5 * kk
    char <- max(c(newr0 - r0, newr1 - r1, newr2 - r2))
    r0 <- newr0
    r1 <- newr1
    r2 <- newr2
    lljilu[itr,5] <- Q(e, y, r0, r1, r2, ECI, hatW, lambda, beta)
    if (itr > minitr) {
      if (chabeta < 0.1 && chalambda < 0.01 && char < 0.01 && chaECI < 0.01) {
        break
      }
      if ((lljilu[itr] - lljilu[itr - 1]) < 0.01 ) {
        break
      }
    }
    resultprint <- paste("The Q function value is ", lljilu[itr], sep = "")
    print(resultprint)
  }
  ## Output result
  colnames(ECI) <- typename
  rownames(ECI) <- typename
  typeresult <- ECI
  cellresult <- data.frame(cellname[cellpair[, 1]], cellname[cellpair[, 2]], e)
  colnames(cellresult) <- c("Cell 1", "Cell 2", "e")
  lljilu <- lljilu[1:itr,]
  result <- list(typeresult, cellresult, lambda, beta, c(r0, r1, r2),lljilu)
  return(result)
}
