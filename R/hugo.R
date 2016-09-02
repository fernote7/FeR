#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @examples
#' cat_function()



#require(RQuantLib)
#require(plot3D)
#require(dplyr)

euler <- function(X0, mu, sigma, Dt, t, T, N, plt=F){ 
    
    X0 = rep(X0,N)
    linhas = (T-t)/Dt + 1
    w = matrix(nrow=linhas, ncol = N)
    n = seq(from=Dt, to=T, by=Dt)
    w[1,] = X0
    
    
    for (j in 1:N){ 
        for (i in seq_along(n)){
        w[i+1,j] = w[i,j] + mu *w[i,j] *Dt + sigma * w[i,j] * sqrt(Dt)*rnorm(1,0,1)
        }
    }
    
    
    if (plt == T){        
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 6){
            for (i in seq_along(w[1,])){
                points(w[,i], col='red', pch=20)
            }
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1)
    }
    lista = list('solution' = as.data.frame(w), 'mean' = mean(w[linhas,]))
    return(invisible(lista))
}



milstein <- function(X0, mu, sigma, Dt, t, T, N, plt=F){ 
    
    X0 = rep(X0,N)
    linhas = (T-t)/Dt + 1
    w = matrix(nrow=linhas, ncol = N)
    n = seq(from=Dt, to=T, by=Dt)
    w[1,] = X0
    
    
    for (j in 1:N){ 
        for (i in seq_along(n)){
            w[i+1,j] = w[i,j] + mu *w[i,j] *Dt + sigma * w[i,j] * sqrt(Dt)*rnorm(1,0,1) +
                        (1/2) * sigma * (Dt * rnorm(1,0,1)^2 - Dt)
        }
    }
    
    
    if (plt == T){    
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 5){
            for (i in seq_along(w[1,])){ 
                points(w[,i], col='red', pch=20)
            }
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1) 
    }
    lista = list('solution' = as.data.frame(w), 'mean' = mean(w[linhas,]))
    return(invisible(lista))
}



vec_euler <- function(X0, mu, sigma, Dt, t, T, N, plt=FALSE){ 
    
    X0 = rep(X0,N)                          # criando o vetor de valores iniciais
    linhas = floor((T-t)/Dt + 1)            # número de linhas da matriz de resultados
    w = matrix(nrow=linhas, ncol = N)       # matriz com N colunas, cada simulação é uma coluna
    n = seq(from=Dt, to=T, by=Dt)           # passos a serem dados
    w[1,] = X0                              # primeira linha da matriz recebe o valor inicial X0
    l2 = linhas-1                           # linhas da matriz de números aleatórios
    bw = sqrt(Dt)*rnorm(N*l2,0,1)           # números aleatórios criados
    bw = matrix(bw, nrow = l2, ncol = N)    # coloca em forma matricial
    
    
    for (i in seq_along(n)){
        w[i+1,] = w[i,] + mu *w[i,] *Dt + sigma * w[i,] * bw[i,]
    }
    
    
    if (plt == TRUE){
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 6){
            for (i in seq_along(w[1,])){
                points(w[,i], col='red', pch=20)
            }
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1)
    }
    lista = list('realizations' = as.data.frame(w), 'mean' = mean(w[linhas,]),
                 'solution' = w[linhas,])
    return(invisible(lista))
}




langevin <- function(X0, mu, sigma, Dt, t, T, N, plt=F){ 
    
    X0 = rep(X0,N)
    linhas = (T-t)/Dt + 1
    w = matrix(nrow=linhas, ncol = N)
    n = seq(from=Dt, to=T, by=Dt)
    w[1,] = X0
    
    
    for (j in 1:N){ 
        for (i in seq_along(n)){
            w[i+1,j] = w[i,j] - mu *w[i,j] *Dt + sigma * sqrt(Dt)*rnorm(1,0,1)
        }
    }
    
    
    if (plt == T){    
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 5){
            for (i in seq_along(w[1,])){ 
                points(w[,i], col='red', pch=20)
            }
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1) 
    }    
    lista = list('solution' = as.data.frame(w))
    return(invisible(lista))
}







sorg <- function(X0, mu, sigma, Dt, t, T, N, plt=F){
    
    
    X0 = rep(X0,N)                          # criando o vetor de valores iniciais
    linhas = floor((T-t)/Dt + 1)            # número de linhas da matriz de resultados
    w = matrix(nrow=linhas, ncol = N)       # matriz com N colunas, cada simulação é uma coluna
    n = seq(from=Dt, to=T, by=Dt)           # passos a serem dados
    w[1,] = X0                              # primeira linha da matriz recebe o valor inicial X0
    l2 = linhas-1                           # linhas da matriz de números aleatórios
    bw = sqrt(Dt)*rnorm(N*l2,0,1)           # números aleatórios criados
    bw = matrix(bw, nrow = l2, ncol = N)    # coloca em forma matricial
    a <- function(x) (mu)
    b <- function(x) (sigma)
        
    
    for (i in seq_along(n)){
        w[i+1,] = w[i,] + mu *w[i,] *Dt + sigma * w[i,] * bw[i,] +
                  (1/2) * (b(w[i,]+b(w[i,])*sqrt(Dt)) - b(w[i,])) * 
                  (Dt * rnorm(1,0,1)^2 - Dt)/sqrt(Dt)
    }
    
    
    if (plt == T){
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 6){
            for (i in seq_along(w[1,])){
                points(w[,i], col='red', pch=20)
            }
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1)
    }
    lista = list('solution' = as.data.frame(w), 'mean' = mean(w[linhas,]))
    return(invisible(lista))
    
    
    
    
}




sorg2 <- function(X0, mu1, sigma, Dt, t, T, N, plt=F){
    
    X0 = rep(X0,N)                          # criando o vetor de valores iniciais
    linhas = floor((T-t)/Dt + 1)            # número de linhas da matriz de resultados
    w = matrix(nrow=linhas, ncol = N)       # matriz com N colunas, cada simulação é uma coluna
    n = seq(from=Dt, to=T, by=Dt)           # passos a serem dados
    w[1,] = X0                              # primeira linha da matriz recebe o valor inicial X0
    l2 = linhas-1                           # linhas da matriz de números aleatórios
    bw = sqrt(Dt)*rnorm(N*l2,0,1)           # números aleatórios criados
    bw = matrix(bw, nrow = l2, ncol = N)    # coloca em forma matricial
    mu <- function(mu1) (eval(mu1))
    b <- function(x) (sigma)
    
    
    for (i in seq_along(n)){
        w[i+1,] = w[i,] + mu *w[i,] *Dt + sigma * w[i,] * bw[i,] +
            (1/2) * mu
    }
    
    
    if (plt == T){
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 6){
            for (i in seq_along(w[1,])){
                points(w[,i], col='red', pch=20)
            }
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1)
    }
    lista = list('solution' = as.data.frame(w), 'mean' = mean(w[linhas,]))
    return(invisible(lista))
    
    
    
    
}



vec_euler_call <- function(X0, mu, sigma, Dt = 0.1, t, T, N = 100, K = 0, plt=FALSE){ 
    
    X0 = rep(X0,N)                          # criando o vetor de valores iniciais
    linhas = floor((T-t)/Dt + 1)            # número de linhas da matriz de resultados
    w = matrix(nrow=linhas, ncol = N)       # matriz com N colunas, cada simulação é uma coluna
    n = seq(from=Dt, to=T, by=Dt)           # passos a serem dados
    w[1,] = X0                              # primeira linha da matriz recebe o valor inicial X0
    l2 = linhas-1                           # linhas da matriz de números aleatórios
    bw = sqrt(Dt)*rnorm(N*l2,0,1)           # números aleatórios criados
    bw = matrix(bw, nrow = l2, ncol = N)    # coloca em forma matricial

    
    for (i in seq_along(n)){
        w[i+1,] = w[i,] + mu *w[i,] *Dt + sigma * w[i,] * bw[i,]
    }
    
    Result = w[linhas,] - K
    Result[Result <= 0] = 0
    call = exp(-mu*T)*mean(Result)
    
    if (plt == TRUE){
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 6){
            for (i in seq_along(w[1,])){
                points(w[,i], col='red', pch=20)
            }
        }
        if (K != 0){
            abline(h=K, col='red', lwd = 3)
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1)
    }
    lista = list('realizations' = as.data.frame(w), 'mean' = mean(w[linhas,]),
                 'solution' = w[linhas,], 'eurocall' = call)
    return(invisible(lista))
}

