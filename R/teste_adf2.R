#' A Cat Function
#' 
#' adf2 function
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @examples
#' cat_function()


teste_adf2 = function(serie, lag_max){
    
    # Carregando pacotes necessários
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(urca)
    
    # Etapa 1: modelo com tendência - Delta Y_t = a0 + a1*t + g Y_t-1 + somatorio(beta_i Delta Y_t-i) + e_t
    adf = ur.df(serie, type = "trend", lags = lag_max, selectlags = "AIC")
    coef_adf = adf@testreg$coefficients
    
    # H0: g = 0?
    # Rejeitamos H0 e concluímos que a série é estacionária se estatística de teste (tau3) for menor que 
    # o valor crítico (cval) ao nível de significância alpha (no nosso caso escolhemos alpha = 5%).
    # Caso tau3 > cval(5%) não rejeitamos H0 e concluímos que a série NÃO é estacionária.
    if(adf@teststat[1] < adf@cval[1,2]){
        resultado_adf = "série estacionária"
    }else{
        resultado_adf = "série não-estacionária"
        
        # Concluído que a série não é estacionária, devemos verificar se o coeficiente de tendência 
        # determinística é significativa. O teste tem como hipótese nula H0: a1 = 0, dado que g = 0, que é 
        # diferente de testar H0: a1 = g = 0 (teste F). Optamos por realizar somente o PRIMEIRO teste.
        # Rejeitamos H0 se phi3 > cval (5% nesse caso) e concluímos que existe tendência determinística.
        if(coef_adf[which(rownames(coef_adf) == "tt"), "Pr(>|t|)"] < 0.05){ #Não
            
            # Como concluímos pela presença de tendência determinística e que a série não é estacionária, devemos
            # realizar um novo teste para o coeficiente g usando os valores críticos da distribuição normal. 
            # O teste é unilateral com H1: g < 0, então cval(5%) = -1.645.
            if(adf@teststat[1] < -1.645){ #Não
                resultado_adf = "série estacionária"
            }else{resultado_adf = "série não-estacionária"}
        }else{#Sim
            
            # O teste concluiu pela ausência de tendência determinística e, portanto, conduzimos o teste ADF
            # somente com drift.
            # Etapa 2: modelo com drift - Delta Y_t = a0 + g Y_t-1 + somatorio(beta_i Delta Y_t-i) + e_t
            adf = ur.df(serie, type = "drift", lags = lag_max, selectlags = "AIC")
            coef_adf = adf@testreg$coefficients
            
            # Testaremos H0: g = 0, ou seja, se a série tem raiz unitária.
            # Rejeitaremos H0 se tau < cval(5%) e concluiremos que a série é estacionária. Caso contrário 
            # (tau > cval(5%)) não rejeitaremos a hipótese nula e testaremos a significância de a0. 
            if(adf@teststat[1] < adf@cval[1,2]){
                resultado_adf = "série estacionária"
            }else{
                
                # Não rejeitamos H0, ou seja, a série não é estacionária. Testaremos então se o coeficiente a0
                # é significativo conjuntamente com o fato da série ser não-estacionária, ou seja, testaremos
                # a hipótese nula H0: a0 = 0 dado que g = 0 (teste t). Alguns autores sugerem testar H0: a0 = g = 0,
                # o que significa usar o teste F. Como os teste podem discordar em alguns casos optamos por realizar 
                # somente o primeiro teste (H0: a0 = 0 dado que g = 0).
                if(coef_adf[which(rownames(coef_adf) == "(Intercept)"), "Pr(>|t|)"] < 0.05){ #Não
                    
                    # Como concluímos pela presença do drift e que a série não é estacionária, devemos
                    # realizar um novo teste para o coeficiente g usando os valores críticos da distribuição normal. 
                    # O teste é unilateral com H1: g < 0, então cval(5%) = -1.645.
                    if(adf@teststat[1] < -1.645){ #Não
                        resultado_adf = "série estacionária"
                    }else{resultado_adf = "série não-estacionária"}
                    
                }else{ #Sim
                    
                    # O teste concluiu pela ausência de drift e, portanto, conduzimos o teste ADF sem drift ou 
                    # tendência determinística.
                    # Etapa 3: modelo sem drift ou tendência - Delta Y_t = g Y_t-1 + somatorio(beta_i Delta Y_t-i) + e_t
                    adf = ur.df(serie, type = "none", lags = lag_max, selectlags = "AIC")
                    coef_adf = adf@testreg$coefficients
                    
                    if(adf@teststat[1] < adf@cval[1,2]){
                        resultado_adf = "série estacionária"
                    }else{resultado_adf = "série não-estacionária"}
                }
            }
        }
    }
    
    if(resultado_adf == "série estacionária"){aux = 0}else{aux = 1}
    return(aux)
}