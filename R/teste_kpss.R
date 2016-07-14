#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

teste_kpss = function(serie){
    
    # Carregando pacotes necessários
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(urca, zoo)
    
    
    # 'tend_det' é uma variável auxiliar que diz se a tendência determinística foi
    # significativa ou não
    tend_det = NA
    
    # Testando se a série tem tendência determinística ao nível de significância de 5%
    trend = seq(1, length(serie), by = 1)
    lm = lm(serie ~ trend)
    sumario = summary(lm)
    pvalor = sumario$coefficients[which(rownames(sumario$coefficients) == "trend"),"Pr(>|t|)"]
    if(pvalor < 0.05){
        
        # Rejeitamos H0: beta1 = 0 no modelo Y = beta0 + beta1*trend se 'pvalor_trend' for menor que 5%.
        # Nesse caso, devemos trabalhar com a covariável sem tendência determinística.
        serie = ts(lm$residuals, 
                   start = c(year(as.Date(serie)[1]), month(as.Date(serie)[1])), 
                   freq = 12)
        tend_det = 1
    }else{tend_det = 0}
    # OBS: Não lidaremos com quebras estruturais nem relações não-lineares.
    
    # Testada a tendência, conduziremos agora o teste KPSS para verificar se a série é estacionária. 
    # kpss = ur.kpss(serie, type = "mu", use.lag = lags)
    kpss = ur.kpss(serie, type = "mu", lags = "long")
    
    # Resultado do teste KPSS
    if(kpss@teststat < kpss@cval[2]){
        # Não rejeito H0: a série é estacionária
        resultado_kpss = "série estacionária"
    }else{resultado_kpss = "série não-estacionária"}
    
    if(resultado_kpss == "série estacionária"){
        aux = list(estac = 0, tend_det = tend_det)
    }else{aux = list(estac = 1, tend_det = tend_det)}
    
    return(aux)
}