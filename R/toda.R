#' ADF Function
#'
#' This function allows you to run an ADF test through several contiguous columns of a dataframe using urca package.
#' @param x Dataframe or tibble.
#' @param m Series maximum integration order.
#' @param ic Information criteria, from 1 to 4 (AIC, HQ, SC, FPE). Defaults to 1 (AIC).
#' @param plt Plot? Defaults to FALSE.
#' @keywords Toda-Yamamoto
#' @export
#' @import vars
#' @examples
#' toda(x, m = 2, plt = F, ic = 1)
#' @author Fernando Teixeira
#' 

# x = granger_ts; m = 2; plt = T; ic = 1
# toda(x = granger_ts, m = 2, plt = T, ic = 1)

toda <- function(x, m, plt = F, ic = 1){
    
    vetor = c()
    
    for ( i in 1:length(x[1,])){
        
        k = i 
        
        
        for (j in k:length(x[1,])){
            
            
            if(names(x[i]) != names(x[j])){
                
                granger_teste = ts(cbind(x[,i], x[,j]))
                colnames(granger_teste)[1] = colnames(x[i])
                colnames(granger_teste)[2] = colnames(x[j])
                aic = VARselect(granger_teste, lag.max = 12, 
                                type = "const")$selection[ic]
                
                a1=VAR(granger_teste, p= aic, type = "const")
                
                erros = serial.test(a1,lags.bg = aic, type = "BG")
                
                
                for (n in (aic):11){
                    if(erros$serial$p.value < 0.05){
                        a1 = VAR(granger_teste, p= (n + 1), type = "const")
                        erros = serial.test(a1,lags.bg = (n + 1), type = "BG")
                    } 
                }
                
                
                
                aic = a1$p
                

                exog = c()
                
                for (n in 1:m){
                    a=lag(granger_teste[,1],(aic+n))
                    b=lag(granger_teste[,2],(aic+n))
                    exog = cbind(exog, a, b)
                }
                
                var.ex <- VAR(granger_teste, p= aic,type = "const", 
                              exogen = exog)
                
                
                
                causa1 = causality(var.ex,
                                   cause = colnames(granger_teste)[1])$Granger 
                causa2 = causality(var.ex,
                                   cause = colnames(granger_teste)[2])$Granger
                
                
                
                
                if (causa1$p.value < 0.1 & causa2$p.value < 0.1){
                    
                    c1 = strsplit(causa1$method, split=" ")
                    vetor = c(vetor, paste(c1[[1]][4], "<->", c1[[1]][8]))
                    
                    if (plt == T){
                        stab=stability(a1)
                        plot(stab)
                    }
                    
                }
                
                
            }    
        }
    }
    lista = list(causalidade = vetor)
    return(invisible(lista))
}