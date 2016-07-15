teste_dfgls = function(serie, modelo = c("trend", "constant"), lag_max){
    
    # Carregando pacotes necessários
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(urca)
    
    # Teste DF-GLS para verificar se a série é estacionária. 
    dfgls = ur.ers(serie, type = "DF-GLS", model = modelo, lag.max = lag_max)
    
    # Resultado do teste KPSS
    if(dfgls@teststat < dfgls@cval[2]){
        # Rejeito H0: a série não é estacionária ao nível de significância de 5%
        resultado_dfgls = "série estacionária"
    }else{resultado_dfgls = "série não-estacionária"}
    
    if(resultado_dfgls == "série estacionária"){aux = 0}else{aux = 1}
    return(aux)
    
}