teste_kpss2 = function(serie, modelo = c("tau", "mu")){
    

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
