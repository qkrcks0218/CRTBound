logistic <- function(x){
    xT <- -50*as.numeric(x<=-50)+x*as.numeric(abs(x)<50)+50*as.numeric(x>=50)
    return( exp(xT)/(1+exp(xT)) )
}

solution <- function(LF,NNT,h,c,thr){
    # LF     : learner f_{\theta}(X)
    # NNT    : total number of NT/AT/COs
    # h,c    : parameters of the surrogate indicator function
    # thr    : numerical error tolerance
    window <- 0.1
    a <- min(LF) ; b <- max(LF)
    qqt <- NULL
    propt <- NULL
    minneg <- 100
    indic <- 0
    INDICATOR <- 0
    if(NNT==length(LF)){
        return(a)
    } else if (NNT==0){
        return(b)
    } else {
        while( (sum( abs(minneg)<thr)==0) | indic <= 1  ){
            qqn <- c(seq(a,b,window),b)
            propn <- rep(0,length(qqn))
            for(ii in 1:length(qqn)){
                x <- LF-qqn[ii]
                coef <- (1-2*c)/(2*c*h)
                A.EE1 <- as.numeric(coef*(x-h)>50)*50 + as.numeric(coef*(x-h)<=50 & coef*(x-h)>-50)*coef*(x-h) +
                    as.numeric(coef*(x-h) <= (-50) )*(-50)
                B.EE1 <- as.numeric(coef*(x+h)>50)*50 + as.numeric(coef*(x+h)<=50 & coef*(x+h)>-50)*coef*(x+h) +
                    as.numeric(coef*(x+h) <= (-50) )*(-50)
                EE1 <- as.numeric(x>=h)*(1-c*exp(-A.EE1)) + as.numeric(x<h&x>=-h)*(coef*(x+h)*c+c) + as.numeric(x< (-h))*(c*exp(B.EE1))
                propn[ii] <- sum(EE1) - NNT
            }
            a <- max(qqn[propn>0])
            b <- min(qqn[propn<0])
            qqt <- c(qqt,qqn)
            propt <- c(propt,propn)
            window <- window/10
            minneg <- max(propt[propt<=0])
            if((sum( abs(minneg)<thr)>0)){
                INDICATOR <- 1
            }
            indic <- indic + INDICATOR

        }
        A <- cbind(qqt,propt)
        A <- A[A[,2]<=0,]
        A <- A[abs(A[,2])<thr,]
        if(is.null(dim(A))){
            return(A[1])
        } else {
            return(median(A[,1]))
        }
    }

}

Data.Reform <- function(Y,Z,A,C,X,seed=1){
    Data <- list()
    Data$Y <- Y
    Data$Z <- Z
    Data$A <- A
    Data$C <- C
    Data$X <- X
    if( sum(apply((X-1)^2,2,sum)==0)==0 ){
        X <- cbind(1,X)
        colnames(X) <- c("Const",colnames(X))
    } else {
        pos.const <- which(apply((X-1)^2,2,sum)==0)
        if( pos.const>1 ){
            X <- X[,c(pos.const,(1:dim(X)[2])[-pos.const])]
        }
    }

    Data$N <- length(Data$Y)

    Data$J <- length(table(Data$C))
    Data$Class <- unique(Data$C)

    Data$n <- rep(0,Data$N)
    for(jj in 1:length(Data$N)){
        Data$n[jj] <- sum(Data$C==Data$C[jj])
    }
    Data$nc <- rep(0,Data$J)
    for(jj in 1:Data$J){
        Data$nc[jj] <- sum(Data$C==Data$Class[jj])
    }

    Data$Zc <- as.numeric(unique(cbind(Data$C,Data$Z))[,2])
    Data$m <- sum(Data$Zc)

    set.seed(seed)

    Data$error <- cbind( 2*runif(Data$N)*(10^(-10)) - 10^(-10),
                         2*runif(Data$N)*(10^(-10)) - 10^(-10),
                         2*runif(Data$N)*(10^(-10)) - 10^(-10) )
    return(Data)
}


Simulation.Reform <- function(Y0,Y1,Z,A0,A1,C,X,seed=1){
    Data <- list()
    Data$Y0 <- Y0
    Data$Y1 <- Y1
    Data$Y <- Y1*Z + Y0*(1-Z)
    Data$Z <- Z
    Data$A0 <- A0
    Data$A1 <- A1
    Data$A <- A1*Z + A0*(1-Z)
    Data$C <- C
    Data$X <- X
    if( sum(apply((X-1)^2,2,sum)==0)==0 ){
        X <- cbind(1,X)
        colnames(X) <- c("Const",colnames(X))
    } else {
        pos.const <- which(apply((X-1)^2,2,sum)==0)
        if( pos.const>1 ){
            X <- X[,c(pos.const,(1:dim(X)[2])[-pos.const])]
        }
    }

    Data$N <- length(Data$Y)

    Data$J <- length(table(Data$C))
    Data$Class <- unique(Data$C)

    Data$n <- rep(0,Data$N)
    for(jj in 1:length(Data$N)){
        Data$n[jj] <- sum(Data$C==Data$C[jj])
    }
    Data$nc <- rep(0,Data$J)
    for(jj in 1:Data$J){
        Data$nc[jj] <- sum(Data$C==Data$Class[jj])
    }

    Data$Zc <- as.numeric(unique(cbind(Data$C,Data$Z))[,2])
    Data$m <- sum(Data$Zc)

    Data$CompSt <- rep("",Data$N)
    Data$CompSt[(A1==0 & A0==0)] <- "NT"
    Data$CompSt[(A1==1 & A0==1)] <- "AT"
    Data$CompSt[(A1==1 & A0==0)] <- "CO"

    set.seed(seed)

    Data$error <- cbind( 2*runif(Data$N)*(10^(-10)) - 10^(-10),
                         2*runif(Data$N)*(10^(-10)) - 10^(-10),
                         2*runif(Data$N)*(10^(-10)) - 10^(-10) )
    return(Data)
}

ITT <- function(Data,Input.Type="Data"){
    # result : reformed data
    Y <- Data$Y
    Z <- Data$Z
    C <- Data$C
    A <- Data$A
    Class <- Data$Class
    J <- Data$J
    N <- Data$N
    n <- Data$n
    nc <- Data$nc
    Zc <- Data$Zc
    m <- Data$m
    X <- Data$X

    if(Input.Type!="Data"){
        Y0 <- Data$Y0
        Y1 <- Data$Y1
    }

    bY.type1 <- rep(0,J)
    for(jj in 1:J){
        if(Zc[jj]==1){
            bY.type1[jj] <- sum(Y[C==Class[jj]])/sum(nc[Zc==1])
        } else {
            bY.type1[jj] <- sum(Y[C==Class[jj]])/sum(nc[Zc==0])
        }
    }

    hat.type1.t1 <- sum((bY.type1*Zc))
    hat.type1.t0 <- sum(bY.type1*(1-Zc))
    R1 <- rep(0,J)
    for(jj in 1:J){
        if(Zc[jj]==1){
            R1[jj] <- ( sum(Y[C==Class[jj]]) - hat.type1.t1*nc[jj] )*J/N
        } else {
            R1[jj] <- ( sum(Y[C==Class[jj]]) - hat.type1.t0*nc[jj] )*J/N
        }

    }
    cons.var.type1 <- var(R1[Zc==1])/m + var(R1[Zc==0])/(J-m)

    OUT1 <- matrix( c(hat.type1.t1 - hat.type1.t0,
                      sqrt(cons.var.type1),
                      (hat.type1.t1-hat.type1.t0)/sqrt(cons.var.type1)),1,3)
    colnames(OUT1) <- c("ITT","SE","z-statistic")
    rownames(OUT1) <- ""


    if(Input.Type!="Data"){

        bY1.type1 <- rep(0,J)
        bY0.type1 <- rep(0,J)
        for(jj in 1:J){
            bY1.type1[jj] <- sum(Y1[C==Class[jj]])*(J/N)
            bY0.type1[jj] <- sum(Y0[C==Class[jj]])*(J/N)
        }

        type1.t1 <- sum(bY1.type1)/J
        type1.t0 <- sum(bY0.type1)/J

        trueR1 <- (bY1.type1 - type1.t1*nc*J/N)
        trueR0 <- (bY0.type1 - type1.t0*nc*J/N)

        true.var.type1 <- var(trueR1)/m + var(trueR0)/(J-m) - var(trueR1-trueR0)/J

        OUT1 <- matrix(c(OUT1,c(type1.t1- type1.t0, sqrt(true.var.type1))),1,5)

        colnames(OUT1) <- c("ITT","SE","z-statistic","True.ITT","True.SE")
        rownames(OUT1) <- ""

    }


    return(OUT1)

}

HTE <- function(Data,Xvar=NULL,constant=TRUE,Input.Type="Data"){
    # result : reformed data
    # Xvar   : column indices corresponding to beta
    Y <- Data$Y
    Z <- Data$Z
    C <- Data$C
    A <- Data$A
    Class <- Data$Class
    J <- Data$J
    N <- Data$N
    n <- Data$n
    nc <- Data$nc
    Zc <- Data$Zc
    m <- Data$m
    if(is.null(Xvar)){
        X <- Data$X
    } else {
        X <- Data$X[,Xvar]
    }

    non.trivial <- which( apply(X,2,var) > 0 )
    X.temp <- X[,non.trivial]
    if(constant==TRUE){
        X <- cbind(1,X.temp)
        colnames(X) <- c("Const",colnames(X.temp))
    } else {
        X <- X.temp
    }

    if(Input.Type!="Data"){
        Y0 <- Data$Y0
        Y1 <- Data$Y1
    }



    NCmat <- matrix(rep(n,dim(X)[2]),dim(X)[1],dim(X)[2])
    N1 <- sum(nc[Zc==1])
    N0 <- sum(nc[Zc==0])

    beta11 <- solve( t(X*Z)%*%(X*Z) ) %*% ( t(X)%*%(Z*Y) )
    beta10 <- solve( t(X*(1-Z))%*%(X*(1-Z)) ) %*% ( t(X)%*%((1-Z)*Y) )
    beta.type1 <- beta11 - beta10

    Sxx1 <- t(X*Z)%*%(X*Z) / N1
    Sxx0 <- t(X*(1-Z))%*%(X*(1-Z)) / N0
    Sxx <- t(X)%*%X/N

    Residual1 <- Y - X%*%beta11
    Residual0 <- Y - X%*%beta10
    RX1 <- matrix(0,J,dim(X)[2])
    RX0 <- matrix(0,J,dim(X)[2])
    for(ii in 1:J){
        RX1[ii,] <- t(X[C==Class[ii],])%*%Residual1[C==Class[ii]]*(J/N)
        RX0[ii,] <- t(X[C==Class[ii],])%*%Residual0[C==Class[ii]]*(J/N)
    }
    RX1 <- RX1[Zc==1,] ; RX0 <- RX0[Zc==0,]
    beta.var.type11 <- solve(Sxx1)%*%(var(RX1)/m)%*%solve(Sxx1) + solve(Sxx0)%*%(var(RX0)/(J-m))%*%solve(Sxx0)

    if(Input.Type!="Data"){

        truebeta11 <- solve( t(X)%*%X ) %*% ( t(X)%*%Y1 )
        truebeta10 <- solve( t(X)%*%X ) %*% ( t(X)%*%Y0 )
        true.beta.type1 <- truebeta11 - truebeta10
        trueResidual1 <- Y1 - X%*%truebeta11
        trueResidual0 <- Y0 - X%*%truebeta10
        trueRX1 <- matrix(0,J,dim(X)[2])
        trueRX0 <- matrix(0,J,dim(X)[2])
        for(ii in 1:J){
            trueRX1[ii,] <- t(X[C==Class[ii],])%*%trueResidual1[C==Class[ii]]*(J/N)
            trueRX0[ii,] <- t(X[C==Class[ii],])%*%trueResidual0[C==Class[ii]]*(J/N)
        }
        true.beta.var.type1 <- solve(Sxx)%*%( var(trueRX1)/m + var(trueRX0)/(J-m) - var(trueRX1-trueRX0)/J )%*%solve(Sxx)

    }

    Result <- list()

    Result$Estimate <- rbind(as.numeric(beta.type1),
                             sqrt(diag(beta.var.type11)),
                             as.numeric(beta.type1)/sqrt(diag(beta.var.type11)))
    colnames(Result$Estimate) <- colnames(X)
    rownames(Result$Estimate) <- c("Estimate","SE","z-statistic")

    if(constant==TRUE){
        Result$NonConst.Chi.Sq.Statistic <- c( t(beta.type1[-1])%*%solve(beta.var.type11[-1,-1])%*%beta.type1[-1],
                                               dim(X)[2]-1)
    } else {
        Result$NonConst.Chi.Sq.Statistic <- c( t(beta.type1)%*%solve(beta.var.type11)%*%beta.type1,
                                               dim(X)[2])
    }
    Result$NonConst.Chi.Sq.Statistic <- data.frame(matrix(Result$NonConst.Chi.Sq.Statistic,1,2))
    colnames(Result$NonConst.Chi.Sq.Statistic) <- c("Statistic","df")
    rownames(Result$NonConst.Chi.Sq.Statistic) <- ""


    if(Input.Type!="Data"){
        Est.temp <- rbind(as.numeric(true.beta.type1),
                          sqrt(diag(true.beta.var.type1)))
        Result$Estimate <- rbind(Result$Estimate,Est.temp)
        colnames(Result$Estimate) <- colnames(X)
        rownames(Result$Estimate) <- c("Estimate","SE","z-statistic","True.Estimate","True.SE")

    }

    return(Result)
}

#' @export
SharpBound <- function(Data , paraC=NULL, method="Linear", CIcalc=FALSE, SSsize=1000, level=0.95, seed=1, Input.Type="Data"){
    # result : reformed data
    # paraC  : column index of X that are used in classifier
    # method : linear (=1) or logistic (=2)
    # CIcalc : none (=0), bootstrap (="Boot")
    # Bsize  : number of bootstrap iterations

    defaultW <- getOption("warn")
    options(warn = -1)

    alpha <- 1-level

    Result <- list()

    maxY <- max(Data$Y)
    minY <- min(Data$Y)
    Y <- (Data$Y-minY)/(maxY-minY)
    Z <- Data$Z
    C <- Data$C
    A <- Data$A
    Class <- Data$Class
    J <- Data$J
    N <- Data$N
    n <- Data$n
    nc <- Data$nc
    Zc <- Data$Zc
    m <- Data$m

    if(is.null(paraC)){
        X <- Data$X
    } else {
        X <- Data$X[,paraC]
    }
    non.trivial <- which( apply(X,2,var) > 0 )
    X.temp <- X[,non.trivial]
    X <- cbind(1,X.temp)
    colnames(X) <- c("Const",colnames(X.temp))

    error <- Data$error

    ### population level

    if(Input.Type!="Data"){

        Y0 <- Data$Y0
        Y1 <- Data$Y1
        CS <- Data$CompSt


        NT <- as.numeric(CS=="NT")
        AT <- as.numeric(CS=="AT")
        CO <- as.numeric(CS=="CO")

        TNT <- sum(CS=="NT")
        TAT <- sum(CS=="AT")
        TCO <- sum(CS=="CO")

        S1 <- sum(Y1)
        S0 <- sum(Y0)

        T1NT <- sum(Y1*NT)
        T1AT <- sum(Y1*AT)
        T1CO <- sum(Y1*CO)

        T0NT <- sum(Y0*NT)
        T0AT <- sum(Y0*AT)
        T0CO <- sum(Y0*CO)

        ## NT classifier

        if(method=="Linear"){
            theta <- as.numeric(lm(NT~X[,-1])$coefficients)
            theta[is.na(theta)] <- 0
            NT.lf <- X%*%theta
        } else if (method=="Logistic"){
            tempD <- cbind(NT,1,X[,-1])
            GLM <- glmnet::glmnet(tempD[,-1],tempD[,1], family="binomial",alpha=0,lambda=0.0001)
            theta <- as.numeric(c(GLM$a0,GLM$beta[-1]))
            theta[is.na(theta)] <- 0
            NT.lf <- predict(GLM,tempD[,-1],type="response")
            NT.lf.link <- predict(GLM,tempD[,-1],type="link")
        }
        NT.newlf <- NT.lf+error[,1]

        qL <- sort(NT.newlf)[N-TNT]
        qU <- sort(NT.newlf)[N-TNT+1]

        h <- (qU-qL)/3/10
        c <- 1/max(log(TNT),log(N-TNT))/10

        thr <- 0.000001
        q <- try(solution(NT.newlf,TNT,h,c,thr),silent=TRUE)
        NT.C <- as.numeric(NT.newlf>=q)

        ## AT classifier

        if(method=="Linear"){
            theta <- as.numeric(lm(AT~X[,-1])$coefficients)
            theta[is.na(theta)] <- 0
            AT.lf <- X%*%theta
        } else if (method=="Logistic"){
            tempD <- cbind(AT,1,X[,-1])
            GLM <- glmnet::glmnet(tempD[,-1],tempD[,1], family="binomial",alpha=0,lambda=0.0001)
            theta <- as.numeric(c(GLM$a0,GLM$beta[-1]))
            theta[is.na(theta)] <- 0
            AT.lf <- predict(GLM,tempD[,-1],type="response")
            AT.lf.link <- predict(GLM,tempD[,-1],type="link")
        }
        AT.newlf <- AT.lf+error[,2]

        qL <- sort(AT.newlf)[N-TAT]
        qU <- sort(AT.newlf)[N-TAT+1]

        h <- (qU-qL)/3/10
        c <- 1/max(log(TAT),log(N-TAT))/10

        thr <- 0.000001
        q <- try(solution(AT.newlf,TAT,h,c,thr),silent=TRUE)
        AT.C <- as.numeric(AT.newlf>=q)

        ## CO classifier

        if(method=="Linear"){
            CO.lf <- -(TNT/N)*NT.lf-(TAT/N)*AT.lf
        } else if (method=="Logistic"){
            CO.lf <- logistic(-(TNT/N)*NT.lf.link-(TAT/N)*AT.lf.link)
        }

        CO.newlf <- CO.lf+error[,3]

        qL <- sort(CO.newlf)[N-TCO]
        qU <- sort(CO.newlf)[N-TCO+1]

        h <- (qU-qL)/3/10
        c <- 1/max(log(TCO),log(N-TCO))/10

        thr <- 0.000001
        q <- try(solution(CO.newlf,TCO,h,c,thr),silent=TRUE)
        CO.C <- as.numeric(CO.newlf>=q)

        T1CNT <- sum(Y1*NT.C)
        T1CAT <- sum(Y1*AT.C)
        T1CCO <- sum(Y1*CO.C)

        T0CNT <- sum(Y0*NT.C)
        T0CAT <- sum(Y0*AT.C)
        T0CCO <- sum(Y0*CO.C)

        RNT <- sum(NT.C*(1-NT))
        RAT <- sum(AT.C*(1-AT))
        RCO <- sum(CO.C*(1-CO))

        TruePara <- data.frame(matrix(c(TNT,TAT,TCO,S1,S0,T1NT,T0AT,T1CNT,T1CAT,T1CCO,T0CNT,T0CAT,T0CCO,RNT,RAT,RCO),1,16))
        colnames(TruePara) <- c("TNT","TAT","TCO","S1","S0","T1NT","T0AT","T1CNT","T1CAT","T1CCO","T0CNT","T0CAT","T0CCO","RNT","RAT","RCO")
    }




    ################ Estimation

    p1 <- m/J ; p0 <- (J-m)/J

    bY <- bNT <- bAT <- rep(0,J)
    bT1NT <- bT0AT <- rep(0,J)
    for(jj in 1:J){
        if(Zc[jj]==1){
            bY[jj] <- sum(Y[C==Class[jj]])/sum(nc[Zc==1])
            bNT[jj] <- sum((1-A)[C==Class[jj]])/sum(nc[Zc==1])
            bT1NT[jj] <- sum((Y*(1-A))[C==Class[jj]])/sum(nc[Zc==1])
        } else {
            bY[jj] <- sum(Y[C==Class[jj]])/sum(nc[Zc==0])
            bAT[jj] <- sum(A[C==Class[jj]])/sum(nc[Zc==0])
            bT0AT[jj] <- sum((Y*A)[C==Class[jj]])/sum(nc[Zc==0])
        }
    }

    hat.S1 <- sum(bY*Zc)*N
    hat.S0 <- sum(bY*(1-Zc))*N

    hat.TNT <- sum(bNT*Zc)*N
    hat.TAT <- sum(bAT*(1-Zc))*N
    hat.TCO <- N - hat.TNT - hat.TAT

    hat.T1NT <- sum(bT1NT*Zc)*N
    hat.T0AT <- sum(bT0AT*(1-Zc))*N


    ## NT classifier

    NTest <- (1-A)*Z

    if(method=="Linear"){
        theta <- as.numeric(lm(NTest[Z==1]~X[Z==1,-1])$coefficients)
        theta[is.na(theta)] <- 0
        hat.NT.lf <- X%*%theta
    } else if (method=="Logistic"){
        tempD <- cbind(NTest,1,X[,-1])
        GLM <- glmnet::glmnet(tempD[Z==1,-1],tempD[Z==1,1], family="binomial",alpha=0,lambda=0.0001)
        theta <- as.numeric(c(GLM$a0,GLM$beta[-1]))
        theta[is.na(theta)] <- 0
        hat.NT.lf <- predict(GLM,tempD[,-1],type="response")
        hat.NT.lf.link <- predict(GLM,tempD[,-1],type="link")
    }
    hat.NT.newlf <- hat.NT.lf+error[,1]

    qL <- sort(hat.NT.newlf)[sum(Z)-sum(NTest)]
    qU <- sort(hat.NT.newlf)[sum(Z)-sum(NTest)+1]

    h <- (qU-qL)/3/10
    c <- 1/max(log(hat.TNT),log(N-hat.TNT))/10

    thr <- 0.000001
    q <- try(solution(hat.NT.newlf[Z==1],sum(NTest),h,c,thr),silent=TRUE)
    if(class(q)=="try-error"){
        q <- (qL+qU)/2
    }
    hat.NT.C <- as.numeric(hat.NT.newlf>=q)

    ## AT classifier

    ATest <- A*(1-Z)

    if(method=="Linear"){
        theta <- as.numeric(lm(ATest[Z==0]~X[Z==0,-1])$coefficients)
        theta[is.na(theta)] <- 0
        hat.AT.lf <- X%*%theta
    } else if (method=="Logistic"){
        tempD <- cbind(ATest,1,X[,-1])
        GLM <- glmnet::glmnet(tempD[Z==0,-1],tempD[Z==0,1], family="binomial",alpha=0,lambda=0.0001)
        theta <- as.numeric(c(GLM$a0,GLM$beta[-1]))
        theta[is.na(theta)] <- 0
        hat.AT.lf <- predict(GLM,tempD[,-1],type="response")
        hat.AT.lf.link <- predict(GLM,tempD[,-1],type="link")
    }

    hat.AT.newlf <- hat.AT.lf+error[,2]

    qL <- sort(hat.AT.newlf)[sum(1-Z)-sum(ATest)]
    qU <- sort(hat.AT.newlf)[sum(1-Z)-sum(ATest)+1]

    h <- (qU-qL)/3/10
    c <- 1/max(log(hat.TAT),log(N-hat.TAT))/10

    thr <- 0.000001
    q <- try(solution(hat.AT.newlf[Z==0],sum(ATest),h,c,thr),silent=TRUE)
    if(class(q)=="try-error"){
        q <- (qL+qU)/2
    }
    hat.AT.C <- as.numeric(hat.AT.newlf>=q)


    ## CO classifier

    if(method=="Linear"){
        hat.CO.lf <- -(hat.TNT/N)*hat.NT.lf-(hat.TAT/N)*hat.AT.lf
    } else if (method=="Logistic"){
        hat.CO.lf <- logistic(-(hat.TNT/N)*hat.NT.lf.link-(hat.TAT/N)*hat.AT.lf.link)
    }
    hat.CO.newlf <- hat.CO.lf+error[,3]

    qL <- sort(hat.CO.newlf)[N-round(hat.TCO)]
    qU <- sort(hat.CO.newlf)[N-round(hat.TCO)+1]

    h <- (qU-qL)/3/10
    c <- 1/max(log(hat.TCO),log(N-hat.TCO))/10

    thr <- 0.000001
    q <- try(solution(hat.CO.newlf,hat.TCO,h,c,thr),silent=TRUE)
    if(class(q)=="try-error"){
        q <- (qL+qU)/2
    }
    hat.CO.C <- as.numeric(hat.CO.newlf>=q)


    bT1CNT <- bT1CAT <- bT1CCO <- rep(0,J)
    bT0CNT <- bT0CAT <- bT0CCO <- rep(0,J)
    bRNT <- bRAT <- bRCO <- rep(0,J)

    nNT <- nAT <- nCO <- rep(0,J)
    for(jj in 1:J){
        nNT[jj] <- sum(hat.NT.C[C==Class[jj]])
        nAT[jj] <- sum(hat.AT.C[C==Class[jj]])
        nCO[jj] <- sum(hat.CO.C[C==Class[jj]])
    }

    for(jj in 1:J){
        if(Zc[jj]==1){
            if(sum(nNT[Zc==1])>0){
                bT1CNT[jj] <- sum((Y*hat.NT.C)[C==Class[jj]])/sum(nNT[Zc==1])
                bRNT[jj] <- sum((hat.NT.C*(1-NTest))[C==Class[jj]])/sum(nNT[Zc==1])
            } else {
                bT1CNT[jj] <- 0
                bRNT[jj] <- 0
            }
            if(sum(nAT[Zc==1])>0){
                bT1CAT[jj] <- sum((Y*hat.AT.C)[C==Class[jj]])/sum(nAT[Zc==1])
            } else {
                bT1CAT[jj] <- 0
            }
            if(sum(nCO[Zc==1])>0){
                bT1CCO[jj] <- sum((Y*hat.CO.C)[C==Class[jj]])/sum(nCO[Zc==1])
                bRCO[jj] <- sum((hat.CO.C*NTest)[C==Class[jj]])/sum(nCO[Zc==1])
            } else {
                bT1CCO[jj] <- 0
                bRCO[jj] <- 0
            }
        } else {
            if(sum(nNT[Zc==0])>0){
                bT0CNT[jj] <- sum((Y*hat.NT.C)[C==Class[jj]])/sum(nNT[Zc==0])
            } else {
                bT0CNT[jj] <- 0
            }
            if(sum(nAT[Zc==0])>0){
                bT0CAT[jj] <- sum((Y*hat.AT.C)[C==Class[jj]])/sum(nAT[Zc==0])
                bRAT[jj] <- sum((hat.AT.C*(1-ATest))[C==Class[jj]])/sum(nAT[Zc==0])
            } else {
                bT0CAT[jj] <- 0
                bRAT[jj] <- 0
            }
            if(sum(nCO[Zc==0])>0){
                bT0CCO[jj] <- sum((Y*hat.CO.C)[C==Class[jj]])/sum(nCO[Zc==0])
                bRCO[jj] <- sum((hat.CO.C*ATest)[C==Class[jj]])/sum(nCO[Zc==0])
            } else {
                bT0CCO[jj] <- 0
                bRCO[jj] <- 0
            }


        }
    }


    hat.T1CNT <- sum(bT1CNT*Zc)*hat.TNT
    hat.T1CAT <- sum(bT1CAT*Zc)*hat.TAT
    hat.T1CCO <- sum(bT1CCO*Zc)*hat.TCO

    hat.T0CNT <- sum(bT0CNT*(1-Zc))*hat.TNT
    hat.T0CAT <- sum(bT0CAT*(1-Zc))*hat.TAT
    hat.T0CCO <- sum(bT0CCO*(1-Zc))*hat.TCO

    hat.RNT <- sum(bRNT*Zc)*hat.TNT
    hat.RAT <- sum(bRAT*(1-Zc))*hat.TAT
    hat.RCO <- (sum(bRCO*Zc)+sum(bRCO*(1-Zc)))*hat.TCO



    EstPara <- data.frame(matrix(c(hat.TNT,hat.TAT,hat.TCO,hat.S1,hat.S0,hat.T1NT,hat.T0AT,hat.T1CNT,hat.T1CAT,hat.T1CCO,hat.T0CNT,hat.T0CAT,hat.T0CCO,hat.RNT,hat.RAT,hat.RCO),1,16))
    colnames(EstPara) <- c("TNT","TAT","TCO","S1","S0","T1NT","T0AT","T1CNT","T1CAT","T1CCO","T0CNT","T0CAT","T0CCO","RNT","RAT","RCO")

    ####### LP

    ## Unified Restrictions
    ## A0N, A1N, B0N, B1N, C0N, C1N,
    ## A0A, A1A, B0A, B1A, C0A, C1A,
    ## A0C, A1C, B0C, B1C, B0C, B1C

    f.con <- matrix(c( 0,1,0,1,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  1,0,1,0,0,0,  0,0,0,0,0,0,
                       0,1,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,1,0,0,0,1,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,0,0,0,1,
                       1,0,0,0,1,0,  0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  1,0,0,0,1,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  1,0,0,0,1,0,
                       0,1,0,1,0,0,  0,1,0,1,0,0,  0,1,0,1,0,0,
                       1,0,1,0,0,0,  1,0,1,0,0,0,  1,0,1,0,0,0, # end of eq const
                       1,-1,0,0,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,1,-1,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,1,-1, 0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  1,-1,0,0,0,0, 0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,1,-1,0,0, 0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,1,-1, 0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  1,-1,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,1,-1,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,1,-1,
                       0,1,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,1,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,1,0,0,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,1,0,0,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,1,  0,0,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,0,0,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,1,0,0,
                       0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,1)  ,28,18, byrow=T)
    f.dir <- c(rep("==",10),rep("<=",18))


    f.obj.NT <- c(-1,1,-1,1,0,0, rep(0,12))
    f.obj.AT <- c(rep(0,6),-1,1,-1,1,0,0, rep(0,6))
    f.obj.CO <- c(rep(0,12), -1,1,-1,1,0,0)

    if(Input.Type!="Data"){
        f.rhs <- c( T1NT, T0AT,
                    T1CNT, T1CAT, T1CCO,
                    T0CNT, T0CAT, T0CCO,
                    S1, S0,
                    0,0,0,0,0,0,0,0,0,
                    TNT-RNT, RNT, RNT,
                    TAT-RAT, RAT, RAT,
                    TCO-RCO, RCO, RCO )

        LP.solution <- data.frame( matrix( c(lpSolve::lp("min",f.obj.NT,f.con,f.dir,f.rhs)$objval,
                                             lpSolve::lp("max",f.obj.NT,f.con,f.dir,f.rhs)$objval,
                                             lpSolve::lp("min",f.obj.AT,f.con,f.dir,f.rhs)$objval,
                                             lpSolve::lp("max",f.obj.AT,f.con,f.dir,f.rhs)$objval,
                                             lpSolve::lp("min",f.obj.CO,f.con,f.dir,f.rhs)$objval,
                                             lpSolve::lp("max",f.obj.CO,f.con,f.dir,f.rhs)$objval)/
                                               c(TNT,TNT,TAT,TAT,TCO,TCO),1,6) )
        LP.solution[,is.nan(as.numeric(LP.solution))] <- 0
        colnames(LP.solution) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")

    }

    hat.f.rhs <- c( hat.T1NT, hat.T0AT,
                    hat.T1CNT, hat.T1CAT, hat.T1CCO,
                    hat.T0CNT, hat.T0CAT, hat.T0CCO,
                    hat.S1, hat.S0,
                    0,0,0,0,0,0,0,0,0,
                    hat.TNT-hat.RNT, hat.RNT, hat.RNT,
                    hat.TAT-hat.RAT, hat.RAT, hat.RAT,
                    hat.TCO-hat.RCO, hat.RCO, hat.RCO )

    Violation <- lpSolve::lp("max",f.obj.NT,f.con,f.dir,hat.f.rhs)$status

    if(Violation==0){


        hat.LP.solution <- data.frame( matrix( c(lpSolve::lp("min",f.obj.NT,f.con,f.dir,hat.f.rhs)$objval,
                                                 lpSolve::lp("max",f.obj.NT,f.con,f.dir,hat.f.rhs)$objval,
                                                 lpSolve::lp("min",f.obj.AT,f.con,f.dir,hat.f.rhs)$objval,
                                                 lpSolve::lp("max",f.obj.AT,f.con,f.dir,hat.f.rhs)$objval,
                                                 lpSolve::lp("min",f.obj.CO,f.con,f.dir,hat.f.rhs)$objval,
                                                 lpSolve::lp("max",f.obj.CO,f.con,f.dir,hat.f.rhs)$objval)/
                                                   c(hat.TNT,hat.TNT,hat.TAT,hat.TAT,hat.TCO,hat.TCO),1,6) )

    } else {


        f.con.goal <- matrix(c( 0,1,0,1,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0, 1,-1, rep(0,9*2+18),
                                0,0,0,0,0,0,  1,0,1,0,0,0,  0,0,0,0,0,0, rep(0,1*2), 1,-1, rep(0,8*2+18),
                                0,1,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,2*2), 1,-1, rep(0,7*2+18),
                                0,0,0,0,0,0,  0,1,0,0,0,1,  0,0,0,0,0,0, rep(0,3*2), 1,-1, rep(0,6*2+18),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,0,0,0,1, rep(0,4*2), 1,-1, rep(0,5*2+18),
                                1,0,0,0,1,0,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,5*2), 1,-1, rep(0,4*2+18),
                                0,0,0,0,0,0,  1,0,0,0,1,0,  0,0,0,0,0,0, rep(0,6*2), 1,-1, rep(0,3*2+18),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  1,0,0,0,1,0, rep(0,7*2), 1,-1, rep(0,2*2+18),
                                0,1,0,1,0,0,  0,1,0,1,0,0,  0,1,0,1,0,0, rep(0,8*2), 1,-1, rep(0,1*2+18),
                                1,0,1,0,0,0,  1,0,1,0,0,0,  1,0,1,0,0,0, rep(0,9*2), 1,-1, rep(0,18), # end of eq const
                                1,-1,0,0,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2), -1, rep(0,17),
                                0,0,1,-1,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+1), -1, rep(0,16),
                                0,0,0,0,1,-1, 0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+2), -1, rep(0,15),
                                0,0,0,0,0,0,  1,-1,0,0,0,0, 0,0,0,0,0,0, rep(0,10*2+3), -1, rep(0,14),
                                0,0,0,0,0,0,  0,0,1,-1,0,0, 0,0,0,0,0,0, rep(0,10*2+4), -1, rep(0,13),
                                0,0,0,0,0,0,  0,0,0,0,1,-1, 0,0,0,0,0,0, rep(0,10*2+5), -1, rep(0,12),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  1,-1,0,0,0,0, rep(0,10*2+6),-1, rep(0,11),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,1,-1,0,0, rep(0,10*2+7),-1, rep(0,10),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,1,-1, rep(0,10*2+8),-1, rep(0,9),
                                0,1,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+9), -1, rep(0,8),
                                0,0,0,1,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+10), -1, rep(0,7),
                                0,0,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+11), -1, rep(0,6),
                                0,0,0,0,0,0,  0,1,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+12), -1, rep(0,5),
                                0,0,0,0,0,0,  0,0,0,1,0,0,  0,0,0,0,0,0, rep(0,10*2+13), -1, rep(0,4),
                                0,0,0,0,0,0,  0,0,0,0,0,1,  0,0,0,0,0,0, rep(0,10*2+14), -1, rep(0,3),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,0,0,0,0, rep(0,10*2+15), -1, rep(0,2),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,1,0,0, rep(0,10*2+16), -1, rep(0,1),
                                0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,1,  rep(0,10*2+17),-1)  ,28,18+38, byrow=T)
        f.dir <- c(rep("==",10),rep("<=",18))

        hat.f.rhs <- c( hat.T1NT, hat.T0AT,
                        hat.T1CNT, hat.T1CAT, hat.T1CCO,
                        hat.T0CNT, hat.T0CAT, hat.T0CCO,
                        hat.S1, hat.S0,
                        0,0,0,0,0,0,0,0,0,
                        hat.TNT-hat.RNT, hat.RNT, hat.RNT,
                        hat.TAT-hat.RAT, hat.RAT, hat.RAT,
                        hat.TCO-hat.RCO, hat.RCO, hat.RCO )

        f.obj.NT.max <- c(-1,1,-1,1,0,0, rep(0,12), rep(-10000,38))
        f.obj.NT.min <- c(-1,1,-1,1,0,0, rep(0,12), rep(10000,38))
        f.obj.AT.max <- c(rep(0,6),-1,1,-1,1,0,0, rep(0,6), rep(-10000,38))
        f.obj.AT.min <- c(rep(0,6),-1,1,-1,1,0,0, rep(0,6), rep(10000,38))
        f.obj.CO.max <- c(rep(0,12), -1,1,-1,1,0,0, rep(-10000,38))
        f.obj.CO.min <- c(rep(0,12), -1,1,-1,1,0,0, rep(10000,38))

        hat.LP.solution <- data.frame(matrix(c( sum( lpSolve::lp("min",f.obj.NT.min, f.con.goal, f.dir, hat.f.rhs)$solution[1:18]*f.obj.NT ),
                                                sum( lpSolve::lp("max",f.obj.NT.max, f.con.goal, f.dir, hat.f.rhs)$solution[1:18]*f.obj.NT ),
                                                sum( lpSolve::lp("min",f.obj.AT.min, f.con.goal, f.dir, hat.f.rhs)$solution[1:18]*f.obj.AT ),
                                                sum( lpSolve::lp("max",f.obj.AT.max, f.con.goal, f.dir, hat.f.rhs)$solution[1:18]*f.obj.AT ),
                                                sum( lpSolve::lp("min",f.obj.CO.min, f.con.goal, f.dir, hat.f.rhs)$solution[1:18]*f.obj.CO ),
                                                sum( lpSolve::lp("max",f.obj.CO.max, f.con.goal, f.dir, hat.f.rhs)$solution[1:18]*f.obj.CO ))/
                                                 c(hat.TNT,hat.TNT,hat.TAT,hat.TAT,hat.TCO,hat.TCO),1,6))

        hat.LP.solution[1] <- max(hat.LP.solution[1],0)
        hat.LP.solution[3] <- max(hat.LP.solution[3],0)
        hat.LP.solution[5] <- max(hat.LP.solution[5],0)

        hat.LP.solution[2] <- min(hat.LP.solution[2],1)
        hat.LP.solution[4] <- min(hat.LP.solution[4],1)
        hat.LP.solution[6] <- min(hat.LP.solution[6],1)


    }

    hat.LP.solution[,is.nan(as.numeric(hat.LP.solution))] <- 0
    colnames(hat.LP.solution) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")


    if(CIcalc==TRUE){

        tt <- 0

        S.EstPara <- data.frame(matrix(0,SSsize,16))
        colnames(S.EstPara) <- c("TNT","TAT","TCO","S1","S0","T1NT","T0AT","T1CNT","T1CAT","T1CCO","T0CNT","T0CAT","T0CCO","RNT","RAT","RCO")

        S.estimates <- matrix(0,SSsize,6)
        S.error <- rep(0,SSsize)

        for(ss in 1:SSsize){

            S.hat.TCO  <- -1

            while(S.hat.TCO<=0){

                set.seed(seed+ss+tt)

                sA <- 0

                S.method <- method

                selected.1 <- sample(which(Zc==1),m,replace = TRUE)
                selected.0 <- sample(which(Zc==0),J-m,replace = TRUE)

                S.J <- J
                S.m <- length(selected.1)

                Zselected.1 <- which(C==Class[selected.1[1]])
                Zselected.0 <- which(C==Class[selected.0[1]])

                for(ii in 2:length(selected.1)){
                    Zselected.1 <- c(Zselected.1,which(C==Class[selected.1[ii]]))
                }
                for(ii in 2:length(selected.0)){
                    Zselected.0 <- c(Zselected.0,which(C==Class[selected.0[ii]]))
                }

                S.Y <- Y[c(Zselected.1,Zselected.0)]
                S.A <- A[c(Zselected.1,Zselected.0)]
                S.N <- length(c(Zselected.1,Zselected.0))
                S.nc <- nc[c(selected.1,selected.0)]
                S.Jind <- rep(1:S.J,S.nc)
                S.m <- length(selected.1)
                S.Zc <- c(rep(1, length(selected.1)),rep(0, length(selected.0)))
                S.Z <- rep(0,S.N)
                S.Z[ 1:length(Zselected.1) ] <- 1
                S.X <- X[c(Zselected.1,Zselected.0),]

                S.C <- rep(1:S.J,S.nc)
                S.Class <- 1:S.J

                sA <- min(sum(S.A[S.Z==0]),sum((1-S.A)[S.Z==1]))


                if(sA<=1){
                    S.method <- "Linear" # Use Linear
                }


                ################ Estimation

                S.p1 <- S.m/S.J ; S.p0 <- (S.J-S.m)/S.J

                bY <- bNT <- bAT <- rep(0,S.J)
                bT1NT <- bT0AT <- rep(0,S.J)
                for(jj in 1:S.J){
                    if(S.Zc[jj]==1){
                        bY[jj] <- sum(S.Y[S.C==S.Class[jj]])/sum(S.nc[S.Zc==1])
                        bNT[jj] <- sum((1-S.A)[S.C==S.Class[jj]])/sum(S.nc[S.Zc==1])
                        bT1NT[jj] <- sum((S.Y*(1-S.A))[S.C==S.Class[jj]])/sum(S.nc[S.Zc==1])
                    } else {
                        bY[jj] <- sum(S.Y[S.C==S.Class[jj]])/sum(S.nc[S.Zc==0])
                        bAT[jj] <- sum(S.A[S.C==S.Class[jj]])/sum(S.nc[S.Zc==0])
                        bT0AT[jj] <- sum((S.Y*S.A)[S.C==S.Class[jj]])/sum(S.nc[S.Zc==0])
                    }
                }


                S.hat.S1 <- sum(bY*S.Zc)*S.N
                S.hat.S0 <- sum(bY*(1-S.Zc))*S.N

                S.hat.TNT <- sum(bNT*S.Zc)*S.N
                S.hat.TAT <- sum(bAT*(1-S.Zc))*S.N
                S.hat.TCO <- S.N - S.hat.TNT - S.hat.TAT

                S.hat.T1NT <- sum(bT1NT*S.Zc)*S.N
                S.hat.T0AT <- sum(bT0AT*(1-S.Zc))*S.N

                if (S.hat.TCO<=0){
                    tt <- tt+1
                }
            }

            ## NT classifier

            S.NTest <- (1-S.A)*S.Z

            if(S.method=="Linear"){
                theta <- as.numeric(lm(S.NTest[S.Z==1]~S.X[S.Z==1,-1])$coefficients)
                theta[is.na(theta)] <- 0
                S.hat.NT.lf <- S.X%*%theta
            } else if (S.method=="Logistic"){
                tempD <- cbind(S.NTest,1,S.X[,-1])
                GLM <- glmnet::glmnet(tempD[S.Z==1,-1],tempD[S.Z==1,1], family="binomial",alpha=0,lambda=0.0001)
                theta <- as.numeric(c(GLM$a0,GLM$beta[-1]))
                theta[is.na(theta)] <- 0
                S.hat.NT.lf <- predict(GLM,tempD[,-1],type="response")
                S.hat.NT.lf.link <- predict(GLM,tempD[,-1],type="link")
            }
            S.hat.NT.newlf <- S.hat.NT.lf + (2*runif(S.N)*10^(-10) - 10^(-10))

            qL <- sort(S.hat.NT.newlf)[sum(S.Z)-sum(S.NTest)]
            qU <- sort(S.hat.NT.newlf)[sum(S.Z)-sum(S.NTest)+1]

            h <- (qU-qL)/3/10
            c <- 1/max(log(S.hat.TNT),log(S.N-S.hat.TNT))/10

            thr <- 0.000001
            q <- try(solution(S.hat.NT.newlf[S.Z==1],sum(S.NTest),h,c,thr),silent=TRUE)
            if(class(q)=="try-error"){
                q <- (qL+qU)/2
            }
            S.hat.NT.C <- as.numeric(S.hat.NT.newlf>=q)

            ## AT classifier

            S.ATest <- S.A*(1-S.Z)

            if(sA>0){

                if(S.method=="Linear"){
                    theta <- as.numeric(lm(S.ATest[S.Z==0]~S.X[S.Z==0,-1])$coefficients)
                    theta[is.na(theta)] <- 0
                    S.hat.AT.lf <- S.X%*%theta
                } else if (S.method=="Logistic"){
                    tempD <- cbind(S.ATest,1,S.X[,-1])
                    GLM <- glmnet::glmnet(tempD[S.Z==0,-1],tempD[S.Z==0,1], family="binomial",alpha=0,lambda=0.0001)
                    theta <- as.numeric(c(GLM$a0,GLM$beta[-1]))
                    theta[is.na(theta)] <- 0
                    S.hat.AT.lf <- predict(GLM,tempD[,-1],type="response")
                    S.hat.AT.lf.link <- predict(GLM,tempD[,-1],type="link")
                }

                S.hat.AT.newlf <- S.hat.AT.lf + (2*runif(S.N)*10^(-10) - 10^(-10))

                qL <- sort(S.hat.AT.newlf)[sum(1-S.Z)-sum(S.ATest)]
                qU <- sort(S.hat.AT.newlf)[sum(1-S.Z)-sum(S.ATest)+1]

                h <- (qU-qL)/3/10
                c <- 1/max(log(S.hat.TAT),log(S.N-S.hat.TAT))/10

                thr <- 0.000001
                q <- try(solution(S.hat.AT.newlf[S.Z==0],sum(S.ATest),h,c,thr),silent=TRUE)
                if(class(q)=="try-error"){
                    q <- (qL+qU)/2
                }
                S.hat.AT.C <- as.numeric(S.hat.AT.newlf>=q)

            } else {
                S.hat.AT.lf <- rep(0,dim(S.X)[1])
                S.hat.AT.C <- rep(0,dim(S.X)[1])
            }


            ## CO classifier

            if(S.method=="Linear"){
                S.hat.CO.lf <- -(S.hat.TNT/S.N)*S.hat.NT.lf-(S.hat.TAT/S.N)*S.hat.AT.lf
            } else if (S.method=="Logistic"){
                S.hat.CO.lf <- logistic( -(S.hat.TNT/S.N)*S.hat.NT.lf.link-(S.hat.TAT/S.N)*S.hat.AT.lf.link )
            }

            S.hat.CO.newlf <- S.hat.CO.lf + (2*runif(S.N)*10^(-10) - 10^(-10))

            qL <- sort(S.hat.CO.newlf)[S.N-round(S.hat.TCO)]
            qU <- sort(S.hat.CO.newlf)[S.N-round(S.hat.TCO)+1]

            h <- (qU-qL)/3/10
            c <- 1/max(log(S.hat.TCO),log(S.N-S.hat.TCO))/10

            thr <- 0.000001
            q <- try(solution(S.hat.CO.newlf,S.hat.TCO,h,c,thr),silent=TRUE)
            if(class(q)=="try-error"){
                q <- (qL+qU)/2
            }
            S.hat.CO.C <- as.numeric(S.hat.CO.newlf>=q)


            bT1CNT <- bT1CAT <- bT1CCO <- rep(0,S.J)
            bT0CNT <- bT0CAT <- bT0CCO <- rep(0,S.J)
            bRNT <- bRAT <- bRCO <- rep(0,S.J)

            nNT <- nAT <- nCO <- rep(0,J)
            for(jj in 1:S.J){
                nNT[jj] <- sum(S.hat.NT.C[S.C==S.Class[jj]])
                nAT[jj] <- sum(S.hat.AT.C[S.C==S.Class[jj]])
                nCO[jj] <- sum(S.hat.CO.C[S.C==S.Class[jj]])
            }

            for(jj in 1:S.J){
                if(S.Zc[jj]==1){
                    if(sum(nNT[S.Zc==1])>0){
                        bT1CNT[jj] <- sum((S.Y*S.hat.NT.C)[S.C==S.Class[jj]])/sum(nNT[S.Zc==1])
                        bRNT[jj] <- sum((S.hat.NT.C*(1-S.NTest))[S.C==S.Class[jj]])/sum(nNT[S.Zc==1])
                    } else {
                        bT1CNT[jj] <- 0
                        bRNT[jj] <- 0
                    }
                    if(sum(nAT[S.Zc==1])>0){
                        bT1CAT[jj] <- sum((S.Y*S.hat.AT.C)[S.C==S.Class[jj]])/sum(nAT[S.Zc==1])
                    } else {
                        bT1CAT[jj] <- 0
                    }
                    if(sum(nCO[S.Zc==1])>0){
                        bT1CCO[jj] <- sum((S.Y*S.hat.CO.C)[S.C==S.Class[jj]])/sum(nCO[S.Zc==1])
                        bRCO[jj] <- sum((S.hat.CO.C*S.NTest)[S.C==S.Class[jj]])/sum(nCO[S.Zc==1])
                    } else {
                        bT1CCO[jj] <- 0
                        bRCO[jj] <- 0
                    }


                } else {
                    if(sum(nNT[S.Zc==0])>0){
                        bT0CNT[jj] <- sum((S.Y*S.hat.NT.C)[S.C==S.Class[jj]])/sum(nNT[S.Zc==0])
                    } else {
                        bT0CNT[jj] <- 0
                    }
                    if(sum(nAT[S.Zc==0])>0){
                        bT0CAT[jj] <- sum((S.Y*S.hat.AT.C)[S.C==S.Class[jj]])/sum(nAT[S.Zc==0])
                        bRAT[jj] <- sum((S.hat.AT.C*(1-S.ATest))[S.C==S.Class[jj]])/sum(nAT[S.Zc==0])
                    } else {
                        bT0CAT[jj] <- 0
                        bRAT[jj] <- 0
                    }
                    if(sum(nCO[S.Zc==0])>0){
                        bT0CCO[jj] <- sum((S.Y*S.hat.CO.C)[S.C==S.Class[jj]])/sum(nCO[S.Zc==0])
                        bRCO[jj] <- sum((S.hat.CO.C*S.ATest)[S.C==S.Class[jj]])/sum(nCO[S.Zc==0])
                    } else {
                        bT0CCO[jj] <- 0
                        bRCO[jj] <- 0
                    }
                }
            }


            S.hat.T1CNT <- sum(bT1CNT*S.Zc)*S.hat.TNT
            S.hat.T1CAT <- sum(bT1CAT*S.Zc)*S.hat.TAT
            S.hat.T1CCO <- sum(bT1CCO*S.Zc)*S.hat.TCO

            S.hat.T0CNT <- sum(bT0CNT*(1-S.Zc))*S.hat.TNT
            S.hat.T0CAT <- sum(bT0CAT*(1-S.Zc))*S.hat.TAT
            S.hat.T0CCO <- sum(bT0CCO*(1-S.Zc))*S.hat.TCO

            S.hat.RNT <- sum(bRNT*S.Zc)*S.hat.TNT
            S.hat.RAT <- sum(bRAT*(1-S.Zc))*S.hat.TAT
            S.hat.RCO <- (sum(bRCO*S.Zc)+sum(bRCO*(1-S.Zc)))*S.hat.TCO


            S.hat.f.rhs <- c( S.hat.T1NT, S.hat.T0AT,
                              S.hat.T1CNT, S.hat.T1CAT, S.hat.T1CCO,
                              S.hat.T0CNT, S.hat.T0CAT, S.hat.T0CCO,
                              S.hat.S1, S.hat.S0,
                              0,0,0,0,0,0,0,0,0,
                              S.hat.TNT-S.hat.RNT, S.hat.RNT, S.hat.RNT,
                              S.hat.TAT-S.hat.RAT, S.hat.RAT, S.hat.RAT,
                              S.hat.TCO-S.hat.RCO, S.hat.RCO, S.hat.RCO )

            if( lpSolve::lp("min",f.obj.NT,f.con,f.dir,S.hat.f.rhs)$status==0 ){
                S.estimates[ss,] <- c(lpSolve::lp("min",f.obj.NT,f.con,f.dir,S.hat.f.rhs)$objval,
                                      lpSolve::lp("max",f.obj.NT,f.con,f.dir,S.hat.f.rhs)$objval,
                                      lpSolve::lp("min",f.obj.AT,f.con,f.dir,S.hat.f.rhs)$objval,
                                      lpSolve::lp("max",f.obj.AT,f.con,f.dir,S.hat.f.rhs)$objval,
                                      lpSolve::lp("min",f.obj.CO,f.con,f.dir,S.hat.f.rhs)$objval,
                                      lpSolve::lp("max",f.obj.CO,f.con,f.dir,S.hat.f.rhs)$objval)/
                    c(S.hat.TNT,S.hat.TNT,S.hat.TAT,S.hat.TAT,S.hat.TCO,S.hat.TCO)
            } else {
                S.error[ss] <- 1
                f.con.goal <- matrix(c( 0,1,0,1,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0, 1,-1, rep(0,9*2+18),
                                        0,0,0,0,0,0,  1,0,1,0,0,0,  0,0,0,0,0,0, rep(0,1*2), 1,-1, rep(0,8*2+18),
                                        0,1,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,2*2), 1,-1, rep(0,7*2+18),
                                        0,0,0,0,0,0,  0,1,0,0,0,1,  0,0,0,0,0,0, rep(0,3*2), 1,-1, rep(0,6*2+18),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,0,0,0,1, rep(0,4*2), 1,-1, rep(0,5*2+18),
                                        1,0,0,0,1,0,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,5*2), 1,-1, rep(0,4*2+18),
                                        0,0,0,0,0,0,  1,0,0,0,1,0,  0,0,0,0,0,0, rep(0,6*2), 1,-1, rep(0,3*2+18),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  1,0,0,0,1,0, rep(0,7*2), 1,-1, rep(0,2*2+18),
                                        0,1,0,1,0,0,  0,1,0,1,0,0,  0,1,0,1,0,0, rep(0,8*2), 1,-1, rep(0,1*2+18),
                                        1,0,1,0,0,0,  1,0,1,0,0,0,  1,0,1,0,0,0, rep(0,9*2), 1,-1, rep(0,18), # end of eq const
                                        1,-1,0,0,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2), -1, rep(0,17),
                                        0,0,1,-1,0,0, 0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+1), -1, rep(0,16),
                                        0,0,0,0,1,-1, 0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+2), -1, rep(0,15),
                                        0,0,0,0,0,0,  1,-1,0,0,0,0, 0,0,0,0,0,0, rep(0,10*2+3), -1, rep(0,14),
                                        0,0,0,0,0,0,  0,0,1,-1,0,0, 0,0,0,0,0,0, rep(0,10*2+4), -1, rep(0,13),
                                        0,0,0,0,0,0,  0,0,0,0,1,-1, 0,0,0,0,0,0, rep(0,10*2+5), -1, rep(0,12),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  1,-1,0,0,0,0, rep(0,10*2+6),-1, rep(0,11),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,1,-1,0,0, rep(0,10*2+7),-1, rep(0,10),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,1,-1, rep(0,10*2+8),-1, rep(0,9),
                                        0,1,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+9), -1, rep(0,8),
                                        0,0,0,1,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+10), -1, rep(0,7),
                                        0,0,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+11), -1, rep(0,6),
                                        0,0,0,0,0,0,  0,1,0,0,0,0,  0,0,0,0,0,0, rep(0,10*2+12), -1, rep(0,5),
                                        0,0,0,0,0,0,  0,0,0,1,0,0,  0,0,0,0,0,0, rep(0,10*2+13), -1, rep(0,4),
                                        0,0,0,0,0,0,  0,0,0,0,0,1,  0,0,0,0,0,0, rep(0,10*2+14), -1, rep(0,3),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,0,0,0,0, rep(0,10*2+15), -1, rep(0,2),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,1,0,0, rep(0,10*2+16), -1, rep(0,1),
                                        0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,1,  rep(0,10*2+17),-1)  ,28,18+38, byrow=T)
                f.dir <- c(rep("==",10),rep("<=",18))

                f.obj.NT.max <- c(-1,1,-1,1,0,0, rep(0,12), rep(-10000,38))
                f.obj.NT.min <- c(-1,1,-1,1,0,0, rep(0,12), rep(10000,38))
                f.obj.AT.max <- c(rep(0,6),-1,1,-1,1,0,0, rep(0,6), rep(-10000,38))
                f.obj.AT.min <- c(rep(0,6),-1,1,-1,1,0,0, rep(0,6), rep(10000,38))
                f.obj.CO.max <- c(rep(0,12), -1,1,-1,1,0,0, rep(-10000,38))
                f.obj.CO.min <- c(rep(0,12), -1,1,-1,1,0,0, rep(10000,38))

                S.estimates[ss,] <- c( sum( lpSolve::lp("min",f.obj.NT.min, f.con.goal, f.dir, S.hat.f.rhs)$solution[1:18]*f.obj.NT ),
                                       sum( lpSolve::lp("max",f.obj.NT.max, f.con.goal, f.dir, S.hat.f.rhs)$solution[1:18]*f.obj.NT ),
                                       sum( lpSolve::lp("min",f.obj.AT.min, f.con.goal, f.dir, S.hat.f.rhs)$solution[1:18]*f.obj.AT ),
                                       sum( lpSolve::lp("max",f.obj.AT.max, f.con.goal, f.dir, S.hat.f.rhs)$solution[1:18]*f.obj.AT ),
                                       sum( lpSolve::lp("min",f.obj.CO.min, f.con.goal, f.dir, S.hat.f.rhs)$solution[1:18]*f.obj.CO ),
                                       sum( lpSolve::lp("max",f.obj.CO.max, f.con.goal, f.dir, S.hat.f.rhs)$solution[1:18]*f.obj.CO ))/
                    c(S.hat.TNT,S.hat.TNT,S.hat.TAT,S.hat.TAT,S.hat.TCO,S.hat.TCO)

            }

            S.estimates[ss,is.nan(S.estimates[ss,])] <- 0
            S.EstPara[ss,] <- c(S.hat.TNT,S.hat.TAT,S.hat.TCO,S.hat.S1,S.hat.S0,S.hat.T1NT,S.hat.T0AT,S.hat.T1CNT,S.hat.T1CAT,S.hat.T1CCO,S.hat.T0CNT,S.hat.T0CAT,S.hat.T0CCO,S.hat.RNT,S.hat.RAT,S.hat.RCO)
        }


    }

    S.estimates[,1] <- apply(cbind(1,apply(cbind(0,S.estimates[,1]),1,max)),1,min)
    S.estimates[,2] <- apply(cbind(1,apply(cbind(0,S.estimates[,2]),1,max)),1,min)
    S.estimates[,3] <- apply(cbind(1,apply(cbind(0,S.estimates[,3]),1,max)),1,min)
    S.estimates[,4] <- apply(cbind(1,apply(cbind(0,S.estimates[,4]),1,max)),1,min)
    S.estimates[,5] <- apply(cbind(1,apply(cbind(0,S.estimates[,5]),1,max)),1,min)
    S.estimates[,6] <- apply(cbind(1,apply(cbind(0,S.estimates[,6]),1,max)),1,min)


    Result$Bound <- hat.LP.solution*(maxY-minY) + minY
    rownames(Result$Bound) <- "Est.Bound"
    colnames(Result$Bound) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")

    if(Input.Type!="Data"){
        Result$Bound <- rbind(hat.LP.solution*(maxY-minY) + minY,
                              LP.solution*(maxY-minY) + minY)
        rownames(Result$Bound) <- c("Est.Bound","True.Bound")
        colnames(Result$Bound) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")
    }

    # Result$Para <- EstPara
    if(Violation==0){
        Result$Violation <- "No Violation"
    } else {
        Result$Violation <- lpSolve::lp("min",f.obj.NT.min, f.con.goal, f.dir, hat.f.rhs)$solution[-(1:18)]
    }

    if(CIcalc==TRUE){
        # Result$Resample.Para <- S.EstPara*(maxY-minY) + minY
        # Result$Resample.Violation <- S.error
        Result$Resample.Bound <- S.estimates*(maxY-minY) + minY
        Result$BootCIBound <- matrix(0,1,6)
        Result$BootCIBound[1] <- max(quantile(S.estimates[,1]*(maxY-minY) + minY,alpha/2),minY)
        Result$BootCIBound[3] <- max(quantile(S.estimates[,3]*(maxY-minY) + minY,alpha/2),minY)
        Result$BootCIBound[5] <- max(quantile(S.estimates[,5]*(maxY-minY) + minY,alpha/2),minY)
        Result$BootCIBound[2] <- min(quantile(S.estimates[,2]*(maxY-minY) + minY,1-alpha/2),maxY)
        Result$BootCIBound[4] <- min(quantile(S.estimates[,4]*(maxY-minY) + minY,1-alpha/2),maxY)
        Result$BootCIBound[6] <- min(quantile(S.estimates[,6]*(maxY-minY) + minY,1-alpha/2),maxY)
        Result$BootCIBound <- data.frame(Result$BootCIBound)
        colnames(Result$Resample.Bound) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")
        colnames(Result$BootCIBound) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")

    }

    options(warn = defaultW)

    return(Result)

}

LongHudgens <- function(Data , paraC=NULL, CIcalc=FALSE, SSsize=1000, level=0.95, seed=1, Input.Type="Data"){
    # result : reformed data
    # paraC  : column index of X that are used in sharpening bounds
    # CIcalc : none (=0), bootstrap (="Boot")
    # Bsize  : number of bootstrap iterations

    defaultW <- getOption("warn")
    options(warn = -1)

    alpha <- 1-level

    Result <- list()

    maxY <- max(Data$Y)
    minY <- min(Data$Y)
    Y <- (Data$Y-minY)/(maxY-minY)
    Z <- Data$Z
    C <- Data$C
    A <- Data$A
    Class <- Data$Class
    J <- Data$J
    N <- Data$N
    n <- Data$n
    nc <- Data$nc
    Zc <- Data$Zc
    m <- Data$m

    X <- Data$X

    if(!is.null(paraC)){

        for(jj in 1:length(paraC)){
            if(length( unique(X[,paraC[jj]]) )!=2){
                print("Chosen Xs must be binary. Use different paraC arguments.")
                return()
            }
        }

    }



    error <- Data$error

    ### population level

    if(Input.Type!="Data"){

        Y0 <- Data$Y0
        Y1 <- Data$Y1
        CS <- Data$CompSt

        NT <- as.numeric(CS=="NT")
        AT <- as.numeric(CS=="AT")
        CO <- as.numeric(CS=="CO")

        TNT <- sum(CS=="NT")
        TAT <- sum(CS=="AT")
        TCO <- sum(CS=="CO")

        S1 <- sum(Y1)
        S0 <- sum(Y0)

        T1NT <- sum(Y1*NT)
        T1AT <- sum(Y1*AT)
        T1CO <- sum(Y1*CO)

        T0NT <- sum(Y0*NT)
        T0AT <- sum(Y0*AT)
        T0CO <- sum(Y0*CO)

        pi0 <- (T0CO+T0NT)/(TCO+TNT)
        gamma <- (TNT)/(TCO+TNT)

        lambda1 <- (T1AT+T1CO)/(TAT+TCO)
        delta <- (TAT)/(TAT+TCO)

        NT.B <- c(max(0,T1NT/TNT - pi0/gamma),min(T1NT/TNT,T1NT/TNT + (-pi0+(1-gamma))/gamma))
        AT.B <- c(max(0,(lambda1-1+delta)/delta-T0AT/TAT),min(1-T0AT/TAT,lambda1/delta-T0AT/TAT))
        CO.B <- c(max(0,(lambda1-delta)/(1-delta)-pi0/(1-gamma)),min(1,lambda1/(1-delta)+(-pi0+gamma)/(1-gamma)))

        if(!is.null(paraC)){

            ########## Binary variable


            Strata.Var <- as.matrix( X[,paraC] )
            Strata.Var <- as.numeric( apply(matrix(2^((dim(Strata.Var)[2]-1):0),dim(Strata.Var)[1],dim(Strata.Var)[2],byrow=T)*as.matrix( X[,paraC] ),1,sum) )+1
            Strata.Var.lv <- sort(unique(Strata.Var))

            NT.S <- AT.S <- CO.S <- list()
            T1NT.S <- T1AT.S <- T1CO.S <- T0NT.S <- T0AT.S <- T0CO.S <- TNT.S <- TAT.S <- TCO.S <- rep(0,length(Strata.Var.lv))
            pi0.S <- gamma.S <- lambda1.S <- delta.S <- rep(0,length(Strata.Var.lv))
            NT.B.S <- AT.B.S <- CO.B.S <- matrix(0,length(Strata.Var.lv),2)

            for(tt in 1:length(Strata.Var.lv)){
                St <- Strata.Var.lv[tt]
                NT.S[[tt]] <- NT*as.numeric(Strata.Var==St)
                AT.S[[tt]] <- AT*as.numeric(Strata.Var==St)
                CO.S[[tt]] <- CO*as.numeric(Strata.Var==St)

                T1NT.S[tt] <- sum( Y1*NT.S[[tt]] )
                T1AT.S[tt] <- sum( Y1*AT.S[[tt]] )
                T1CO.S[tt] <- sum( Y1*CO.S[[tt]] )

                T0NT.S[tt] <- sum( Y0*NT.S[[tt]] )
                T0AT.S[tt] <- sum( Y0*AT.S[[tt]] )
                T0CO.S[tt] <- sum( Y0*CO.S[[tt]] )

                TNT.S[tt] <- sum(NT.S[[tt]])
                TAT.S[tt] <- sum(AT.S[[tt]])
                TCO.S[tt] <- sum(CO.S[[tt]])

                pi0.S[tt] <- (T0CO.S[tt]+T0NT.S[tt])/(TCO.S[tt]+TNT.S[tt])
                gamma.S[tt] <- TNT.S[tt]/(TCO.S[tt]+TNT.S[tt])
                lambda1.S[tt] <- (T1AT.S[tt]+T1CO.S[tt])/(TAT.S[tt]+TCO.S[tt])
                delta.S[tt] <- TAT.S[tt]/(TAT.S[tt]+TCO.S[tt])

                NT.B.S[tt,] <- c(max(0,T1NT.S[tt]/TNT.S[tt] - pi0.S[tt]/gamma.S[tt]),
                                 min(T1NT.S[tt]/TNT.S[tt],T1NT.S[tt]/TNT.S[tt] + (-pi0.S[tt]+(1-gamma.S[tt]))/gamma.S[tt]))
                AT.B.S[tt,] <- c(max(0,(lambda1.S[tt]-1+delta.S[tt])/delta.S[tt]-T0AT.S[tt]/TAT.S[tt]),
                                 min(1-T0AT.S[tt]/TAT.S[tt],lambda1.S[tt]/delta.S[tt]-T0AT.S[tt]/TAT.S[tt]))
                CO.B.S[tt,] <- c(max(0,(lambda1.S[tt]-delta.S[tt])/(1-delta.S[tt])-pi0.S[tt]/(1-gamma.S[tt])),
                                 min(1,lambda1.S[tt]/(1-delta.S[tt])+(-pi0.S[tt]+gamma.S[tt])/(1-gamma.S[tt])))

                NT.B.S[tt,is.nan(NT.B.S[tt,])] <- 0
                AT.B.S[tt,is.nan(AT.B.S[tt,])] <- 0
                CO.B.S[tt,is.nan(CO.B.S[tt,])] <- 0

            }

            TNT.S[is.nan(TNT.S)] <- 0
            TAT.S[is.nan(TAT.S)] <- 0
            TCO.S[is.nan(TCO.S)] <- 0


            NT.B.A <- apply(matrix(TNT.S/TNT,length(Strata.Var.lv),2)*NT.B.S,2,sum)
            AT.B.A <- apply(matrix(TAT.S/TAT,length(Strata.Var.lv),2)*AT.B.S,2,sum)
            CO.B.A <- apply(matrix(TCO.S/TCO,length(Strata.Var.lv),2)*CO.B.S,2,sum)
        }

    }


    ################ Estimation

    p1 <- m/J ; p0 <- (J-m)/J

    bY <- bNT <- bAT <- rep(0,J)
    bT1NT <- bT0AT <- rep(0,J)
    for(jj in 1:J){
        if(Zc[jj]==1){
            bY[jj] <- sum(Y[C==Class[jj]])/sum(nc[Zc==1])
            bNT[jj] <- sum((1-A)[C==Class[jj]])/sum(nc[Zc==1])
            bT1NT[jj] <- sum((Y*(1-A))[C==Class[jj]])/sum(nc[Zc==1])
        } else {
            bY[jj] <- sum(Y[C==Class[jj]])/sum(nc[Zc==0])
            bAT[jj] <- sum(A[C==Class[jj]])/sum(nc[Zc==0])
            bT0AT[jj] <- sum((Y*A)[C==Class[jj]])/sum(nc[Zc==0])
        }
    }

    hat.S1 <- sum(bY*Zc)*N
    hat.S0 <- sum(bY*(1-Zc))*N

    hat.TNT <- sum(bNT*Zc)*N
    hat.TAT <- sum(bAT*(1-Zc))*N
    hat.TCO <- N - hat.TNT - hat.TAT

    hat.T1NT <- sum(bT1NT*Zc)*N
    hat.T0AT <- sum(bT0AT*(1-Zc))*N


    hat.pi0 <- (hat.S0 - hat.T0AT)/(hat.TCO+hat.TNT)
    hat.gamma <- (hat.TNT)/(hat.TCO+hat.TNT)

    hat.lambda1 <- (hat.S1 - hat.T1NT)/(hat.TAT+hat.TCO)
    hat.delta <- (hat.TAT)/(hat.TAT+hat.TCO)

    hat.NT.B <- c(max(0,hat.T1NT/hat.TNT - hat.pi0/hat.gamma),min(hat.T1NT/hat.TNT,hat.T1NT/hat.TNT + (-hat.pi0+(1-hat.gamma))/hat.gamma))
    hat.AT.B <- c(max(0,(hat.lambda1-1+hat.delta)/hat.delta-hat.T0AT/hat.TAT),min(1-hat.T0AT/hat.TAT,hat.lambda1/hat.delta-hat.T0AT/hat.TAT))
    hat.CO.B <- c(max(0,(hat.lambda1-hat.delta)/(1-hat.delta)-hat.pi0/(1-hat.gamma)),min(1,hat.lambda1/(1-hat.delta)+(-hat.pi0+hat.gamma)/(1-hat.gamma)))



    ########## Binary variable

    if(!is.null(paraC)){

        Strata.Var <- as.matrix( X[,paraC] )
        Strata.Var <- as.numeric( apply(matrix(2^((dim(Strata.Var)[2]-1):0),dim(Strata.Var)[1],dim(Strata.Var)[2],byrow=T)*as.matrix( X[,paraC] ),1,sum) )+1
        Strata.Var.lv <- sort(unique(Strata.Var))

        VarType <- unique( cbind(Strata.Var,X[,paraC]) )
        VarType <- VarType[order(VarType[,1]),]

        N.S <- N.Z1.S <- N.Z0.S <- hat.S1.S <- hat.S0.S <- hat.T1NT.S <- hat.T1AT.S <- hat.T1CO.S <-
            hat.T0NT.S <- hat.T0AT.S <- hat.T0CO.S <- hat.TNT.S <- hat.TAT.S <- hat.TCO.S <- rep(0,length(Strata.Var.lv))
        hat.pi0.S <- hat.gamma.S <- hat.lambda1.S <- hat.delta.S <- rep(0,length(Strata.Var.lv))
        hat.NT.B.S <- hat.AT.B.S <- hat.CO.B.S <- matrix(0,length(Strata.Var.lv),2)

        for(tt in 1:length(Strata.Var.lv)){
            St <- Strata.Var.lv[tt]
            xx <- as.numeric(Strata.Var==St)
            N.S[tt] <- sum(xx)
            N.Z1.S[tt] <- sum(xx*as.numeric(Z==1))
            N.Z0.S[tt] <- sum(xx*as.numeric(Z==0))

            bY <- bNT <- bAT <- rep(0,J)
            bT1NT <- bT0AT <- rep(0,J)
            for(jj in 1:J){
                if(Zc[jj]==1){
                    if( sum(xx[Z==1])!=0 ){
                        bY[jj] <- sum((Y*xx)[C==Class[jj]])/sum(xx[Z==1])
                        bNT[jj] <- sum(((1-A)*xx)[C==Class[jj]])/sum(xx[Z==1])
                        bT1NT[jj] <- sum((Y*(1-A)*xx)[C==Class[jj]])/sum(xx[Z==1])
                    } else {
                        bY[jj] <- bNT[jj] <- bT1NT[jj] <- 0
                    }
                } else {
                    if( sum(xx[Z==0])!= 0 ){
                        bY[jj] <- sum((Y*xx)[C==Class[jj]])/sum(xx[Z==0])
                        bAT[jj] <- sum((A*xx)[C==Class[jj]])/sum(xx[Z==0])
                        bT0AT[jj] <- sum((Y*A*xx)[C==Class[jj]])/sum(xx[Z==0])
                    } else {
                        bY[jj] <- bAT[jj] <- bT0AT[jj] <- 0
                    }

                }
            }

            hat.S1.S[tt] <- sum(bY*Zc)*N.S[tt]
            hat.S0.S[tt] <- sum(bY*(1-Zc))*N.S[tt]

            hat.TNT.S[tt] <- sum(bNT*Zc)*N.S[tt]
            hat.TAT.S[tt] <- sum(bAT*(1-Zc))*N.S[tt]
            hat.TCO.S[tt] <- max(N.S[tt] - hat.TNT.S[tt] - hat.TAT.S[tt],0)

            hat.T1NT.S[tt] <- sum(bT1NT*Zc)*N.S[tt]
            hat.T0AT.S[tt] <- sum(bT0AT*(1-Zc))*N.S[tt]


            hat.pi0.S[tt] <- max(0,min(1,(hat.S0.S[tt] - hat.T0AT.S[tt])/(hat.TCO.S[tt]+hat.TNT.S[tt])))
            hat.gamma.S[tt] <- max(0,min(1,hat.TNT.S[tt]/(hat.TCO.S[tt]+hat.TNT.S[tt])))
            hat.lambda1.S[tt] <- max(0,min(1,(hat.S1.S[tt] - hat.T1NT.S[tt])/(hat.TAT.S[tt]+hat.TCO.S[tt])))
            hat.delta.S[tt] <- max(0,min(1,hat.TAT.S[tt]/(hat.TAT.S[tt]+hat.TCO.S[tt])))

            hat.NT.B.S[tt,] <- c(max(0,hat.T1NT.S[tt]/hat.TNT.S[tt] - hat.pi0.S[tt]/hat.gamma.S[tt]),
                                 min(hat.T1NT.S[tt]/hat.TNT.S[tt],hat.T1NT.S[tt]/hat.TNT.S[tt] + (-hat.pi0.S[tt]+(1-hat.gamma.S[tt]))/hat.gamma.S[tt]))
            hat.AT.B.S[tt,] <- c(max(0,(hat.lambda1.S[tt]-1+hat.delta.S[tt])/hat.delta.S[tt]-hat.T0AT.S[tt]/hat.TAT.S[tt]),
                                 min(1-hat.T0AT.S[tt]/hat.TAT.S[tt],hat.lambda1.S[tt]/hat.delta.S[tt]-hat.T0AT.S[tt]/hat.TAT.S[tt]))
            hat.CO.B.S[tt,] <- c(max(0,(hat.lambda1.S[tt]-hat.delta.S[tt])/(1-hat.delta.S[tt])-hat.pi0.S[tt]/(1-hat.gamma.S[tt])),
                                 min(1,hat.lambda1.S[tt]/(1-hat.delta.S[tt])+(-hat.pi0.S[tt]+hat.gamma.S[tt])/(1-hat.gamma.S[tt])))

            hat.NT.B.S[tt,is.nan(hat.NT.B.S[tt,])] <- 0
            hat.AT.B.S[tt,is.nan(hat.AT.B.S[tt,])] <- 0
            hat.CO.B.S[tt,is.nan(hat.CO.B.S[tt,])] <- 0
        }
        hat.TNT.S[is.nan(hat.TNT.S)] <- 0
        hat.TAT.S[is.nan(hat.TAT.S)] <- 0
        hat.TCO.S[is.nan(hat.TCO.S)] <- 0

        hat.NT.B.A <- apply(matrix(hat.TNT.S/sum(hat.TNT.S),length(Strata.Var.lv),2)*hat.NT.B.S,2,sum)
        hat.AT.B.A <- apply(matrix(hat.TAT.S/sum(hat.TAT.S),length(Strata.Var.lv),2)*hat.AT.B.S,2,sum)
        hat.CO.B.A <- apply(matrix(hat.TCO.S/sum(hat.TCO.S),length(Strata.Var.lv),2)*hat.CO.B.S,2,sum)



    }


    if(CIcalc==TRUE){

        S.estimates <- matrix(0,SSsize,6)
        S.estimates.A <- matrix(0,SSsize,6)


        for(ss in 1:SSsize){

            set.seed(seed+ss)

            sA <- 0

            while(sA<=1){

                selected.1 <- sample(which(Zc==1),m,replace = TRUE)
                selected.0 <- sample(which(Zc==0),J-m,replace = TRUE)

                S.J <- J
                S.m <- length(selected.1)

                Zselected.1 <- which(C==Class[selected.1[1]])
                Zselected.0 <- which(C==Class[selected.0[1]])

                for(ii in 2:length(selected.1)){
                    Zselected.1 <- c(Zselected.1,which(C==Class[selected.1[ii]]))
                }
                for(ii in 2:length(selected.0)){
                    Zselected.0 <- c(Zselected.0,which(C==Class[selected.0[ii]]))
                }

                S.Y <- Y[c(Zselected.1,Zselected.0)]
                S.A <- A[c(Zselected.1,Zselected.0)]
                S.N <- length(c(Zselected.1,Zselected.0))
                S.nc <- nc[c(selected.1,selected.0)]
                S.Jind <- rep(1:S.J,S.nc)
                S.m <- length(selected.1)
                S.Zc <- c(rep(1, length(selected.1)),rep(0, length(selected.0)))
                S.Z <- rep(0,S.N)
                S.Z[ 1:length(Zselected.1) ] <- 1
                S.X <- X[c(Zselected.1,Zselected.0),]

                S.C <- rep(1:S.J,S.nc)
                S.Class <- 1:S.J

                sA <- min(sum(S.A[S.Z==0]),sum((1-S.A)[S.Z==1]))
            }

            ################ Estimation

            S.p1 <- S.m/S.J ; S.p0 <- (S.J-S.m)/S.J

            bY <- bNT <- bAT <- rep(0,S.J)
            bT1NT <- bT0AT <- rep(0,S.J)
            for(jj in 1:S.J){
                if(S.Zc[jj]==1){
                    bY[jj] <- sum(S.Y[S.C==S.Class[jj]])/sum(S.nc[S.Zc==1])
                    bNT[jj] <- sum((1-S.A)[S.C==S.Class[jj]])/sum(S.nc[S.Zc==1])
                    bT1NT[jj] <- sum((S.Y*(1-S.A))[S.C==S.Class[jj]])/sum(S.nc[S.Zc==1])
                } else {
                    bY[jj] <- sum(S.Y[S.C==S.Class[jj]])/sum(S.nc[S.Zc==0])
                    bAT[jj] <- sum(S.A[S.C==S.Class[jj]])/sum(S.nc[S.Zc==0])
                    bT0AT[jj] <- sum((S.Y*S.A)[S.C==S.Class[jj]])/sum(S.nc[S.Zc==0])
                }
            }


            S.hat.S1 <- sum(bY*S.Zc)*S.N
            S.hat.S0 <- sum(bY*(1-S.Zc))*S.N

            S.hat.TNT <- sum(bNT*S.Zc)*S.N
            S.hat.TAT <- sum(bAT*(1-S.Zc))*S.N
            S.hat.TCO <- max(S.N - S.hat.TNT - S.hat.TAT,0)

            S.hat.T1NT <- sum(bT1NT*S.Zc)*S.N
            S.hat.T0AT <- sum(bT0AT*(1-S.Zc))*S.N


            S.hat.pi0 <- max(0,min(1,(S.hat.S0 - S.hat.T0AT)/(S.hat.TCO+S.hat.TNT)))
            S.hat.gamma <- max(0,min(1,(S.hat.TNT)/(S.hat.TCO+S.hat.TNT)))

            S.hat.lambda1 <- max(0,min(1,(S.hat.S1 - S.hat.T1NT)/(S.hat.TAT+S.hat.TCO)))
            S.hat.delta <- max(0,min(1,(S.hat.TAT)/(S.hat.TAT+S.hat.TCO)))

            S.hat.NT.B <- c(max(0,S.hat.T1NT/S.hat.TNT - S.hat.pi0/S.hat.gamma),
                            min(S.hat.T1NT/S.hat.TNT,S.hat.T1NT/S.hat.TNT + (-S.hat.pi0+(1-S.hat.gamma))/S.hat.gamma))
            S.hat.AT.B <- c(max(0,(S.hat.lambda1-1+S.hat.delta)/S.hat.delta-S.hat.T0AT/S.hat.TAT),
                            min(1-S.hat.T0AT/S.hat.TAT,S.hat.lambda1/S.hat.delta-S.hat.T0AT/S.hat.TAT))
            S.hat.CO.B <- c(max(0,(S.hat.lambda1-S.hat.delta)/(1-S.hat.delta)-S.hat.pi0/(1-S.hat.gamma)),
                            min(1,S.hat.lambda1/(1-S.hat.delta)+(-S.hat.pi0+S.hat.gamma)/(1-S.hat.gamma)))

            S.estimates[ss,] <- c(S.hat.NT.B, S.hat.AT.B, S.hat.CO.B)


            if(!is.null(paraC)){

                ########## Binary variable


                S.Strata.Var <- as.matrix( S.X[,paraC] )
                S.Strata.Var <- as.numeric( apply(matrix(2^((dim(S.Strata.Var)[2]-1):0),dim(S.Strata.Var)[1],dim(S.Strata.Var)[2],byrow=T)*as.matrix( S.X[,paraC] ),1,sum) )+1
                S.Strata.Var.lv <- sort(unique(S.Strata.Var))


                S.N.S <- S.hat.S1.S <- S.hat.S0.S <- S.hat.T1NT.S <- S.hat.T1AT.S <- S.hat.T1CO.S <-
                    S.hat.T0NT.S <- S.hat.T0AT.S <- S.hat.T0CO.S <- S.hat.TNT.S <- S.hat.TAT.S <- S.hat.TCO.S <- rep(0,length(Strata.Var.lv))
                S.hat.pi0.S <- S.hat.gamma.S <- S.hat.lambda1.S <- S.hat.delta.S <- rep(0,length(Strata.Var.lv))
                S.hat.NT.B.S <- S.hat.AT.B.S <- S.hat.CO.B.S <- matrix(0,length(Strata.Var.lv),2)

                for(tt in 1:length(S.Strata.Var.lv)){
                    S.St <- S.Strata.Var.lv[tt]
                    S.xx <- as.numeric(S.Strata.Var==S.St)
                    S.N.S[tt] <- sum(S.xx)

                    S.bY <- S.bNT <- S.bAT <- rep(0,S.J)
                    S.bT1NT <- S.bT0AT <- rep(0,S.J)
                    for(jj in 1:S.J){
                        if(S.Zc[jj]==1){
                            if(sum(S.xx[S.Z==1])!=0){
                                S.bY[jj] <- sum((S.Y*S.xx)[S.C==S.Class[jj]])/sum(S.xx[S.Z==1])
                                S.bNT[jj] <- sum(((1-S.A)*S.xx)[S.C==S.Class[jj]])/sum(S.xx[S.Z==1])
                                S.bT1NT[jj] <- sum((S.Y*(1-S.A)*S.xx)[S.C==S.Class[jj]])/sum(S.xx[S.Z==1])
                            } else {
                                S.bY[jj] <- S.bNT[jj] <- S.bT1NT[jj] <- 0
                            }
                        } else {
                            if(sum(S.xx[S.Z==0])!=0){
                                S.bY[jj] <- sum((S.Y*S.xx)[S.C==S.Class[jj]])/sum(S.xx[S.Z==0])
                                S.bAT[jj] <- sum((S.A*S.xx)[S.C==S.Class[jj]])/sum(S.xx[S.Z==0])
                                S.bT0AT[jj] <- sum((S.Y*S.A*S.xx)[S.C==S.Class[jj]])/sum(S.xx[S.Z==0])
                            } else {
                                S.bY[jj] <- S.bAT[jj] <- S.bT0AT[jj] <- 0
                            }

                        }
                    }

                    S.hat.S1.S[tt] <- sum(S.bY*S.Zc)*S.N.S[tt]
                    S.hat.S0.S[tt] <- sum(S.bY*(1-S.Zc))*S.N.S[tt]

                    S.hat.TNT.S[tt] <- sum(S.bNT*S.Zc)*S.N.S[tt]
                    S.hat.TAT.S[tt] <- sum(S.bAT*(1-S.Zc))*S.N.S[tt]
                    S.hat.TCO.S[tt] <- max( S.N.S[tt] - S.hat.TNT.S[tt] - S.hat.TAT.S[tt], 0)

                    S.hat.T1NT.S[tt] <- sum(S.bT1NT*S.Zc)*S.N.S[tt]
                    S.hat.T0AT.S[tt] <- sum(S.bT0AT*(1-S.Zc))*S.N.S[tt]


                    S.hat.pi0.S[tt] <- max(0,min(1,(S.hat.S0.S[tt] - S.hat.T0AT.S[tt])/(S.hat.TCO.S[tt]+S.hat.TNT.S[tt])))
                    S.hat.gamma.S[tt] <- max(0,min(1,S.hat.TNT.S[tt]/(S.hat.TCO.S[tt]+S.hat.TNT.S[tt])))
                    S.hat.lambda1.S[tt] <- max(0,min(1,(S.hat.S1.S[tt] - S.hat.T1NT.S[tt])/(S.hat.TAT.S[tt]+S.hat.TCO.S[tt])))
                    S.hat.delta.S[tt] <- max(0,min(1,S.hat.TAT.S[tt]/(S.hat.TAT.S[tt]+S.hat.TCO.S[tt])))

                    S.hat.NT.B.S[tt,] <- c(max(0,S.hat.T1NT.S[tt]/S.hat.TNT.S[tt] - S.hat.pi0.S[tt]/S.hat.gamma.S[tt]),
                                           min(S.hat.T1NT.S[tt]/S.hat.TNT.S[tt],S.hat.T1NT.S[tt]/S.hat.TNT.S[tt] + (-S.hat.pi0.S[tt]+(1-S.hat.gamma.S[tt]))/S.hat.gamma.S[tt]))
                    S.hat.AT.B.S[tt,] <- c(max(0,(S.hat.lambda1.S[tt]-1+S.hat.delta.S[tt])/S.hat.delta.S[tt]-S.hat.T0AT.S[tt]/S.hat.TAT.S[tt]),
                                           min(1-S.hat.T0AT.S[tt]/S.hat.TAT.S[tt],S.hat.lambda1.S[tt]/S.hat.delta.S[tt]-S.hat.T0AT.S[tt]/S.hat.TAT.S[tt]))
                    S.hat.CO.B.S[tt,] <- c(max(0,(S.hat.lambda1.S[tt]-S.hat.delta.S[tt])/(1-S.hat.delta.S[tt])-S.hat.pi0.S[tt]/(1-S.hat.gamma.S[tt])),
                                           min(1,S.hat.lambda1.S[tt]/(1-S.hat.delta.S[tt])+(-S.hat.pi0.S[tt]+S.hat.gamma.S[tt])/(1-S.hat.gamma.S[tt])))
                    S.hat.NT.B.S[tt,is.nan(S.hat.NT.B.S[tt,])] <- 0
                    S.hat.AT.B.S[tt,is.nan(S.hat.AT.B.S[tt,])] <- 0
                    S.hat.CO.B.S[tt,is.nan(S.hat.CO.B.S[tt,])] <- 0
                }

                S.hat.TNT.S[is.nan(S.hat.TNT.S)] <- 0
                S.hat.TAT.S[is.nan(S.hat.TAT.S)] <- 0
                S.hat.TCO.S[is.nan(S.hat.TCO.S)] <- 0

                S.hat.NT.B.A <- apply(matrix(S.hat.TNT.S/sum(S.hat.TNT.S),length(Strata.Var.lv),2)*S.hat.NT.B.S,2,sum)
                S.hat.AT.B.A <- apply(matrix(S.hat.TAT.S/sum(S.hat.TAT.S),length(Strata.Var.lv),2)*S.hat.AT.B.S,2,sum)
                S.hat.CO.B.A <- apply(matrix(S.hat.TCO.S/sum(S.hat.TCO.S),length(Strata.Var.lv),2)*S.hat.CO.B.S,2,sum)

                S.hat.NT.B.A[is.nan(S.hat.NT.B.A)] <- 0
                S.hat.AT.B.A[is.nan(S.hat.AT.B.A)] <- 0
                S.hat.CO.B.A[is.nan(S.hat.CO.B.A)] <- 0

                S.estimates.A[ss,] <- c(S.hat.NT.B.A, S.hat.AT.B.A, S.hat.CO.B.A)

            }

        }

    }

    S.estimates <- as.matrix(S.estimates)

    if(!is.null(paraC)){
        S.estimates.A <- as.matrix(S.estimates.A)
    }



    Result$Bound.NoAdj <- data.frame(matrix(c(hat.NT.B,hat.AT.B,hat.CO.B)*(maxY-minY)+minY,1,6))
    rownames(Result$Bound.NoAdj) <- "Est.Bound"
    if(!is.null(paraC)){
        Result$Bound.Adj <- data.frame(matrix(c(hat.NT.B.A,hat.AT.B.A,hat.CO.B.A)*(maxY-minY)+minY,1,6))
        # Result$Para <- cbind( VarType,
        #                       hat.TNT.S,hat.TCO.S,hat.TAT.S,
        #                       N.S,
        #                       N.Z1.S,N.Z0.S)
        rownames(Result$Bound.Adj) <- "Est.Bound"
    }



    if(Input.Type!="Data"){
        Result$Bound.NoAdj <- data.frame(rbind(c(hat.NT.B,hat.AT.B,hat.CO.B)*(maxY-minY)+minY,
                                               c(NT.B,AT.B,CO.B)*(maxY-minY)+minY))
        Result$Bound.Adj <- data.frame(rbind(c(hat.NT.B.A,hat.AT.B.A,hat.CO.B.A)*(maxY-minY)+minY,
                                             c(NT.B.A,AT.B.A,CO.B.A)*(maxY-minY)+minY))
        rownames(Result$Bound.NoAdj) <- rownames(Result$Bound.Adj) <- c("Est.Bound","True.Bound")
    }


    if(CIcalc==TRUE){
        Result$Resample.Bound.NoAdj <- S.estimates*(maxY-minY)+minY
        Result$BootCIBound.NoAdj <- data.frame(matrix(c( quantile(S.estimates[,1]*(maxY-minY)+minY,alpha/2),
                                                         quantile(S.estimates[,2]*(maxY-minY)+minY,1-alpha/2),
                                                         quantile(S.estimates[,3]*(maxY-minY)+minY,alpha/2),
                                                         quantile(S.estimates[,4]*(maxY-minY)+minY,1-alpha/2),
                                                         quantile(S.estimates[,5]*(maxY-minY)+minY,alpha/2),
                                                         quantile(S.estimates[,6]*(maxY-minY)+minY,1-alpha/2) ),1,6))

        if(!is.null(paraC)){

            Result$Resample.Bound.Adj <- S.estimates.A*(maxY-minY)+minY
            Result$BootCIBound.Adj <- data.frame(matrix(c( quantile(S.estimates.A[,1]*(maxY-minY)+minY,alpha/2),
                                                           quantile(S.estimates.A[,2]*(maxY-minY)+minY,1-alpha/2),
                                                           quantile(S.estimates.A[,3]*(maxY-minY)+minY,alpha/2),
                                                           quantile(S.estimates.A[,4]*(maxY-minY)+minY,1-alpha/2),
                                                           quantile(S.estimates.A[,5]*(maxY-minY)+minY,alpha/2),
                                                           quantile(S.estimates.A[,6]*(maxY-minY)+minY,1-alpha/2) ),1,6))

        }

        colnames(Result$Bound.NoAdj) <-
            colnames(Result$Resample.Bound.NoAdj) <-
            colnames(Result$BootCIBound.NoAdj) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")

        if(!is.null(paraC)){
            colnames(Result$Bound.Adj) <-
                colnames(Result$Resample.Bound.Adj) <-
                colnames(Result$BootCIBound.Adj) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")
        }

    }

    options(warn = defaultW)

    return(Result)

}

Bound.Intersect <- function(SharpBound.Object,LongHudgens.Object,level=0.95){

    alpha <- 1-level

    Result <- list()

    if(dim(SharpBound.Object$Bound)[1]==1){

        Est.Bound <- rep(0,6)
        Est.Bound[c(1,3,5)] <- apply(rbind( SharpBound.Object$Bound[1,] , LongHudgens.Object$Bound.Adj[1,] ),2,max)[c(1,3,5)]
        Est.Bound[c(2,4,6)] <- apply(rbind( SharpBound.Object$Bound[1,] , LongHudgens.Object$Bound.Adj[1,] ),2,min)[c(2,4,6)]
        Result$Bound <- data.frame(matrix(Est.Bound,1,6))
        colnames(Result$Bound) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")

    } else {

        Est.Bound <- matrix(0,2,6)
        Est.Bound[1,c(1,3,5)] <- apply(rbind( SharpBound.Object$Bound[1,] , LongHudgens.Object$Bound.Adj[1,] ),2,max)[c(1,3,5)]
        Est.Bound[1,c(2,4,6)] <- apply(rbind( SharpBound.Object$Bound[1,] , LongHudgens.Object$Bound.Adj[1,] ),2,min)[c(2,4,6)]
        Est.Bound[2,c(1,3,5)] <- apply(rbind( SharpBound.Object$Bound[2,] , LongHudgens.Object$Bound.Adj[2,] ),2,max)[c(1,3,5)]
        Est.Bound[2,c(2,4,6)] <- apply(rbind( SharpBound.Object$Bound[2,] , LongHudgens.Object$Bound.Adj[2,] ),2,min)[c(2,4,6)]
        Result$Bound <- data.frame(Est.Bound)
        colnames(Result$Bound) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")
        rownames(Result$Bound) <- c("Est.Bound","True.Bound")

    }



    if(!is.null(SharpBound.Object$BootCIBound)){

        bootCI <- rep(0,6)

        if(SharpBound.Object$Bound[1,1]>=LongHudgens.Object$Bound.Adj[1,1]){ bootCI[1] <- SharpBound.Object$BootCIBound[1] } else { bootCI[1] <- LongHudgens.Object$BootCIBound.Adj[1] }
        if(SharpBound.Object$Bound[1,3]>=LongHudgens.Object$Bound.Adj[1,3]){ bootCI[3] <- SharpBound.Object$BootCIBound[3] } else { bootCI[3] <- LongHudgens.Object$BootCIBound.Adj[3] }
        if(SharpBound.Object$Bound[1,5]>=LongHudgens.Object$Bound.Adj[1,5]){ bootCI[5] <- SharpBound.Object$BootCIBound[5] } else { bootCI[5] <- LongHudgens.Object$BootCIBound.Adj[5] }
        if(SharpBound.Object$Bound[1,2]<=LongHudgens.Object$Bound.Adj[1,2]){ bootCI[2] <- SharpBound.Object$BootCIBound[2] } else { bootCI[2] <- LongHudgens.Object$BootCIBound.Adj[2] }
        if(SharpBound.Object$Bound[1,4]<=LongHudgens.Object$Bound.Adj[1,4]){ bootCI[4] <- SharpBound.Object$BootCIBound[4] } else { bootCI[4] <- LongHudgens.Object$BootCIBound.Adj[4] }
        if(SharpBound.Object$Bound[1,6]<=LongHudgens.Object$Bound.Adj[1,6]){ bootCI[6] <- SharpBound.Object$BootCIBound[6] } else { bootCI[6] <- LongHudgens.Object$BootCIBound.Adj[6] }

        Result$BootCIBound <- data.frame(matrix(as.numeric(bootCI),1,6))
        colnames(Result$BootCIBound) <- c("NT.LB","NT.UB","AT.LB","AT.UB","CO.LB","CO.UB")

    }

    return(Result)
}






