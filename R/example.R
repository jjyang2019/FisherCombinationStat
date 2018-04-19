
example <- function(){
    set.seed(1)
    pval <- matrix(runif(10*6),10,6)
    pheno <- data.frame(y1=rnorm(20),
                        y2=rnorm(20),
                        y3=rnorm(20) > 0,
                        y4=sample(c("M","F"),20,replace=TRUE),
                        y5=ordered(cut(rnorm(20),4)),
                        y6=ordered(cut(rnorm(20),5))
                        )
    Fisher.comb.pval(pval,pheno,
                     kendall.only=FALSE,miss.pval=FALSE)

    rho <- 0.5
    SIGMA <- matrix(c(1,rho,rho,1),2,2)
    z <-MASS::mvrnorm(n=1000,mu=c(0,0),Sigma=SIGMA)
    x <- as.numeric(cut(z[,1],4))
    y <- as.numeric(cut(z[,2],4))


}
