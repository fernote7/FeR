#' kpss2 function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats


teste_kpss2 = function(serie, modelo = c("tau", "mu")){
    
    # Carregando pacotes necessários
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(urca)
    
    # Teste KPSS para verificar se a série é estacionária
    kpss = ur.kpss(serie, type = modelo, lags = "long")
    
    # Resultado do teste KPSS
    if(kpss@teststat < kpss@cval[2]){
        # Não rejeito H0: a série é estacionária ao nível de significância de 5%
        resultado_kpss = "série estacionária"
    }else{resultado_kpss = "série não-estacionária"}
    
    if(resultado_kpss == "série estacionária"){aux = 0}else{aux = 1}
    return(aux)
}
