
bvpsPDE = function(S0 = 100, r = 0.05, sigma = 0.2, T = 0.25, K1 = 90, K2 = 110, Mt = 300, Mx = 300, smin = 0.000001, smax = 200){
  
  xmin <- log(smin)
  xmax <- log(smax)
  t <- rep(0, times = Mt)
  x <- rep(0, times = Mx)
  
  
  delta_t  <- 0.5*sigma^2*T/(Mt-1)
  for(i in 2:Mt){
    t[i] <- (i-1)*delta_t
  }
  
  delta_x  <- (xmax - xmin)/(Mx-1)
  for(i in 1:Mx){
    x[i]<- xmin +(i-1)*delta_x
  }
  
  q <- (2*r)/(sigma^2)
  
  lambda <- (delta_t)/(delta_x^2)
  
  h <- function(x){
    wyplata <- max((K2-x),0) - max((K1-x),0)
    return(wyplata)
  }
  
  y_initial <- function(x){
    a <- exp(0.5*(q-1)*x)*h(exp(x))
    return(a)
  }
  
  y_initial_xmin <- function(t){
    a <- exp(0.5*(q-1)*xmin + 0.25*(q+1)^2*t)*((K2-K1)*exp(-2*r*t/(sigma^2)))
    return(a)
  }

  prices <- matrix(0,nrow = 2, ncol = Mx)
  
  for(i in 1:Mx){
    prices[1,i] <- y_initial(x[i])
  }
  
  y_initial_xmax <- 0
  
  
  print("Proszê wybraæ schemat")
  var = readline(prompt = "1.Implicit Euler 2.Crank-Nicolson) : ");
  var = as.character(var)
  

  if(var == "1"){
    A = diag((2*lambda+1), ncol = (Mx-2) , nrow =(Mx-2))
    
    for(i in 1:(Mx-3)){
      A[i+1,i] <- -1*lambda
      A[i,i+1]<- -1*lambda
    }
    
    for(i in 1:(Mt-1)){
      prices[2,1] <- lambda* y_initial_xmin(t[i+1])
      prices[2,Mx] <- lambda* y_initial_xmax
      prices[2,(2:(Mx-1))] <- as.numeric(solve(A)%*%prices[1,(2:(Mx-1))])
      prices[1,] <- prices[2,]
      
    }
 
    cena <- rep(0, times = Mx)

    
    for(i in 1:Mx){
      cena[i] <- exp(-0.5*(q-1)*x[i]-0.25*(q+1)^2*t[Mt])*prices[2,i]
    }
    S <- exp(x)
    

    delta <- rep(0, times = (Mx-1))
    for(i in 2:Mx){
      delta[i-1] <- (cena[i]-cena[i-1])/(S[i]-S[i-1])
    }
  
    
    ileft <- max(which(S <= S0))
    iright <-ileft +1 
    
    cenaint <- cena[iright] *(S0-S[ileft])/(S[iright]-S[ileft]) + cena[ileft]*(S[iright]-S0)/(S[iright]-S[ileft])
    
    deltaint <- delta[iright] *(S0-S[ileft])/(S[iright]-S[ileft]) + delta[ileft]*(S[iright]-S0)/(S[iright]-S[ileft])
    
    print(paste("Cena opcji (Implicit Euler): ", cenaint))
    print(paste("Delta opcji (Implicit Euler): ", deltaint))
    
    par(mfrow = c(1,2))
    plot(S,cena, type = "l", xlab = "S",
         ylab = "Cena opcji", col = "blue",
         main = "Wartoœæ opcji (Implicit Euler)")
    abline(v = S0)
    
    wektor_s_short = S[-1]
    index = min(which(delta == 0))
    wektor_s_short2 = wektor_s_short[index:length(wektor_s_short)]
    delta2 = delta[index:length(delta)]
    
    plot(wektor_s_short2,delta2, type = "l", xlab = "S",
         ylab = "Delta", col = "blue",
         main = "Delta opcji (Implicit Euler)",
         ylim = c((min(delta2)-1),(max(delta2)+1)))
    abline(v = S0)
    
  }
  
  
  if(var == "2"){
    
    theta <- 0.5
    S <- exp(x)
    F <- diag(2, ncol = (Mx-2) , nrow =(Mx-2))
    
    for(i in 1:(Mx-3)){
      F[i+1,i] = -1
      F[i,i+1] = -1
    }
    I = diag(1, ncol = (Mx-2) , nrow =(Mx-2))
    
    X <- (I+ theta* lambda* F)
    Y <- (I + (1 - theta)*lambda*F)

    for(i in 1:(Mt-1)){
      prices[2,1] <- lambda * theta* y_initial_xmin(t[i+1]) + (1 - theta) *lambda * y_initial_xmin(t[i])
      prices[2,Mx] <- y_initial_xmax
      prices[2,(2:(Mx-1))] <- as.numeric(solve(X)%*%Y%*%prices[1,(2:(Mx-1))])
      prices[1,]<-prices[2,]
      
    }
    
    cena2 <- rep(0, times = Mx)
    
    for(i in 1:Mx){
      cena2[i] <- exp(-0.5*(q-1)*x[i]-0.25*(q+1)^2*t[Mt])*prices[2,i]
    }
    
    
    delta <- rep(0, times = (Mx-1))
    for(i in 2:Mx){
      delta[i-1] <- (cena2[i]-cena2[i-1])/(S[i]-S[i-1])
    }
    

    left <- max(which(S <= S0))
    right <- left +1 
    
    cenaint <- cena2[right] *(S0-S[left])/(S[right]-S[left]) + cena2[left]*(S[right]-S0)/(S[right]-S[left])
    
    deltaint <- delta[right] *(S0-S[left])/(S[right]-S[left]) + delta[left]*(S[right]-S0)/(S[right]-S[left])
    
    print(paste("Cena opcji (Crank-Nicolson): ", cenaint))
    print(paste("Delta opcji (Crank-Nicolson): ", deltaint))
    
    par(mfrow = c(1,2))
    
    plot(S,cena2, type = "l", xlab = "S",
         ylab = "Cena opcji", col = "blue",
         main = "Wartoœæ opcji (Crank-Nicolson)")
    abline(v = S0)
    
    wektor_s_short = S[-1]
    index = min(which(delta == 0))
    wektor_s_short2 = wektor_s_short[index:length(wektor_s_short)]
    delta2 = delta[index:length(delta)]
    
    plot(wektor_s_short2,delta2, type = "l", xlab = "S",
         ylab = "Delta", col = "blue",
         main = "Delta opcji (Crank-Nicolson)",
         ylim = c((min(delta2)-1),(max(delta2)+1))) 
    abline(v = S0)
  
  }
  if((var !="1") & (var !="2")){
    print("Spróbuj ponownie")
  }
  
}
bvpsPDE(S0=120)

data = read.delim("cp4_data.txt", sep ="=", header = FALSE, dec =".")
data = as.vector(data$V2)

bvpsPDE(100,0.05,0.2,0.25,90,110,200,200,0.00001,200)
bvpsPDE(data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10])
bvpsPDE = function(S0 = 100, r = 0.05, sigma = 0.2, T = 0.25, K1 = 90, K2 = 110, Mt = 200, Mx = 200, smin = 0.00001, smax = 200)




