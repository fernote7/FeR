#' Milstein Function
#'
#' This function allows you to simulate Milstein method.
#' @param X0 Initial value.
#' @param mu Drift coeficient.
#' @param sigma Diffusion coefficient.
#' @param Dt Step size.
#' @param t Initial time period.
#' @param T Final time period.
#' @param N Number of simulations.
#' @param plt Plot? Defaults to FALSE.
#' @keywords milstein method
#' @export
#' @examples
#' @author Fernando Teixeira
#' vec_milstein(307.65, 0.75, 0.3, 0.001, 0, 1, 10000)
#' 
#' 


vec_milstein <- function(X0, mu, sigma, Dt, t, T, N, plt=F){ 
    
    X0 = rep(X0,N)
    linhas = (T-t)/Dt + 1
    w = matrix(nrow=linhas, ncol = N)
    n = seq(from=Dt, to=T, by=Dt)
    w[1,] = X0
    l2 = linhas-1                           # linhas da matriz de números aleatórios
    bw = sqrt(Dt)*rnorm(N*l2,0,1)           # números aleatórios criados
    bw = matrix(bw, nrow = l2, ncol = N)    # coloca em forma matricial
    bw2 = (1/2) * sigma * (Dt * rnorm(N*l2,0,1)^2 - Dt)    # números aleatórios criados 2
    bw2 = matrix(bw2, nrow = l2, ncol = N)                  # coloca em forma matricial 2
    
    
    for (i in seq_along(n)){
        w[i+1,] = w[i,] + mu *w[i,] *Dt + sigma * w[i,] * bw[i,] + bw2[i,]
        
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
    lista = list('realizations' = as.data.frame(w), 'mean' = mean(w[linhas,]),
                 'solution' = w[linhas,])
    return(invisible(lista))
}