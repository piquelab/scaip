### Fun-1
nb_llik <- function(theta, x, size){
## negative binomial log likelihood
mu <- size*exp(theta[1])
inv_disp <- exp(theta[2])

res <- x*log(mu/inv_disp)- x*log(1+mu/inv_disp)-
       inv_disp*log(1+mu/inv_disp)+
       lgamma(x+inv_disp)-lgamma(x+1)-lgamma(inv_disp)
llik <- -sum(res)
} 

###
zinb_llik <- function(theta, x, size){
### log likelihood of x as distributed as ZINB
mu <- size*exp(theta[1])
inv_disp <- exp(theta[2])
logodds <- theta[3]
odds <- exp(logodds)

### if reads equal to zero
if ( sum(x<1)>0){
   x0 <- x[x<1]
   mu0 <- mu[x<1]
   res_zero <- log( odds/(1+odds)+1/(1+odds)*exp(x0*log(mu0/inv_disp)-x0*log(1+mu0/inv_disp)-
               inv_disp*log(1+mu0/inv_disp)+
               lgamma(x0+inv_disp)-lgamma(x0+1)-lgamma(inv_disp)) )
}else{
   res_zero <- 0
}
               
### if reads larger than zero
x1 <- x[x>=1]
mu1 <- mu[x>=1]
res_non_zero <- -log(1+odds)+ 
                x1*log(mu1/inv_disp)- x1*log(1+mu1/inv_disp)-
                inv_disp*log(1+mu1/inv_disp)+
                lgamma(x1+inv_disp)-lgamma(x1+1)-lgamma(inv_disp)
                
llik <- -(sum(res_zero)+sum(res_non_zero))/length(x)
}
