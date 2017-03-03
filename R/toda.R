#' ADF Function
#'
#' This function allows you to run an ADF test through several contiguous columns of a dataframe using urca package.
#' @param x Dataframe or tibble.
#' @param m Series maximum integration order.
#' @param ic Information criteria, from 1 to 4 (AIC, HQ, SC, FPE). Defaults to 1 (AIC).
#' @param plt Plot? Defaults to FALSE.
#' @keywords Toda-Yamamoto
#' @export
#' @import vars aod stats
#' @examples
#' toda(x, m = 2, plt = F, ic = 1, pval=0.05)
#' @author Fernando Teixeira
#' 

# x = granger_ts; m = 2; plt = T; ic = 1, pval=0.1
# yama = toda(x = granger_ts, m = 2, plt = T, ic = 1, pval=0.1)
# boundary(stab$stability[[1]])

toda <- function(x, m, plt = F, ic = 1, pval=0.05){
    
    vetor = c()
    
    for ( i in 1:length(x[1,])){
        
        k = i 
        
        
        for (j in k:length(x[1,])){
            
            
            if(names(x[i]) != names(x[j])){
                
                granger_teste = stats::ts(cbind(x[,i], x[,j]))
                colnames(granger_teste)[1] = colnames(x[i])
                colnames(granger_teste)[2] = colnames(x[j])
                aic = vars::VARselect(granger_teste, lag.max = 12, 
                                type = "const")$selection[ic]
                
                a1=vars::VAR(granger_teste, p= aic, type = "both")
                
                erros = vars::serial.test(a1,lags.bg = aic, type = "BG")
                
                
                ### Decidindo o VAR ideal
                
                for (iter in (aic):12){
                    
                    n = iter
                    
                    if(erros$serial$p.value < 0.05){
                        a1 = vars::VAR(granger_teste, p= (n + 1), type = "both")
                        erros = vars::serial.test(a1,lags.bg = (n + 1), type = "BG")
                    } 
                    
                    stab=vars::stability(a1)
                    b1 = strucchange::boundary(stab$stability[[1]])
                    b2 = -strucchange::boundary(stab$stability[[1]])
                    
                    n2 = n + 1
                    aux = ifelse(
                            (stab[[1]][[1]][[1]] > b1 | stab[[1]][[1]][[1]] < b2) |
                            (stab[[1]][[2]][[1]] > b1 | stab[[1]][[2]][[1]] < b2), 
                                 (n <- n[1]), a1 <- vars::VAR(granger_teste, p= (n2), 
                                                        type = "both"))
                    
                    
                }
                
                rm(aux)
                
                ###
                
                aic = a1$p
                term1 = seq(from = 2, to = (aic*2), by=2)
                term2 = seq(from = 1, to = (aic*2 - 1), by=2)
                
                var.ex <- vars::VAR(granger_teste, p= (aic + m), type = "both")
                rm(a1)
                
                #### Teste de causalidade
                
                causa1 = aod::wald.test(b=coef(var.ex$varresult[[1]]), 
                                   Sigma=vcov(var.ex$varresult[[1]]), Terms=term1)
                
                causa2 = aod::wald.test(b=coef(var.ex$varresult[[2]]), 
                                   Sigma=vcov(var.ex$varresult[[2]]), Terms=term2)
                
                
                
                
                ### Apresentação de resultados
                
                if (causa1$result$chi2[3] < pval & 
                    causa2$result$chi2[3] < pval){
                    
                    c1 = strsplit(colnames(granger_teste), split=" ")
                    vetor = c(vetor, paste(c1[[1]][1], "<->", c1[[2]][1]))
                    
                    if (plt == T){
                        stab2=stability(var.ex)
                        plot(stab2)
                    }
                    
                } else if (causa1$result$chi2[3] < pval) {
                    c1 = strsplit(colnames(granger_teste), split = " ")
                    vetor = c(vetor, paste(c1[[1]][1], "<-", c1[[2]][1]))
                    if (plt == T) {
                        stab2 = stability(var.ex)
                        plot(stab2)
                    }
                } else if (causa2$result$chi2[3] < pval) {
                    c1 = strsplit(colnames(granger_teste), split = " ")
                    vetor = c(vetor, paste(c1[[1]][1], "->", c1[[2]][1]))
                    if (plt == T) {
                        stab2 = stability(var.ex)
                        plot(stab2)
                    }
                }
                
                
                
            }    
        }
    }
    lista = list(causalidade = vetor)
    return(invisible(lista))
}