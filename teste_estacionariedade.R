teste_estacionariedade = function(sondagens, financeiras, periodo, variacao){

    ## Carregando pacotes necessários
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(plyr, lubridate)


    ## Caso o argumento variação esteja faltando consideramos que as séries não estão em %
    if(missing("variacao")){variacao = FALSE}


    ## Pré-tratamento das variáveis: disponibilidade, nulos/negativos nas sondagens
    if(!missing("sondagens")){

        # Caso o argumento 'sondagens' seja incluído, tiramos da base as séries não disponíveis
        # para todo o período e séries de sondagem que contenham valores nulos ou negativos
        # (após ajuste sazonal) no período.

        # Excluindo séries de sondagem não disponíveis no período
        sondagens_aux = sondagens[which(sondagens$Data %in% periodo),]
        aux_exc = is.na(sondagens_aux)
        aux2_exc = apply(aux_exc, 2, sum)
        sondagens_naoexc = sondagens_aux[, names(aux2_exc[aux2_exc == 0])]
        colnames(sondagens_naoexc) = c("Data", paste0("S_", colnames(sondagens_naoexc)[-1]))

        # Caso as séries de sondagem estejam em variação, os valores nulos ou negativos fazem
        # sentido e tais séries não devem ser eliminadas. Para isso incluímos o argumento
        # 'variacao = FALSE'.
        if(variacao == FALSE){

            # Excluindo séries de sondagem com valores negativos ou iguais a zero
            aux3_exc = sondagens_naoexc <= 0
            aux4_exc = apply(aux3_exc, 2, sum)
            sondagens_naoexc = sondagens_naoexc[, names(aux4_exc[aux4_exc == 0])]
        }

        if(!missing(financeiras)){

            # Caso o argumento 'financeiras' seja incluído, tiramos da base as séries financeiras
            # não disponíveis para todo o período

            # Excluindo séries financeiras não disponíveis no período
            financeiras_aux = financeiras[which(financeiras$Data %in% periodo),]
            aux_fin = is.na(financeiras_aux)
            aux2_fin = apply(aux_fin, 2, sum)
            financeiras_naoexc = financeiras_aux[, names(aux2_fin[aux2_fin == 0])]
            colnames(financeiras_naoexc) = c("Data", paste0("F_", colnames(financeiras_naoexc)[-1]))

            # Juntando séries de sondagem e financeiras (base final)
            covariaveis = join_all(list(financeiras_naoexc, sondagens_naoexc),
                                   by = 'Data', type = 'full') # juntando "sondagens" e "financeiras" com base na data
            covariaveis = arrange(covariaveis, Data) # ordenando os dados pela data
            covariaveis_ts = ts(covariaveis[,-1], start = c(year(covariaveis[1,1]), month(covariaveis[1,1])), freq = 12)
        }else{

            # Base final caso somente 'sondagens' seja válido
            covariaveis_ts = ts(sondagens_naoexc[,-1], start = c(year(sondagens_naoexc[1,1]),
                                                                 month(sondagens_naoexc[1,1])), freq = 12)
        }
    }else{

        # Caso 'financeiras' seja incluído, mas 'sondagens' não

        # Excluindo séries financeiras não disponíveis no período
        financeiras_aux = financeiras[which(financeiras$Data %in% periodo),]
        aux_fin = is.na(financeiras_aux)
        aux2_fin = apply(aux_fin, 2, sum)
        financeiras_naoexc = financeiras_aux[, names(aux2_fin[aux2_fin == 0])]
        colnames(financeiras_naoexc) = c("Data", paste0("F_", colnames(financeiras_naoexc)[-1]))

        # Base final
        covariaveis_ts = ts(financeiras_naoexc[,-1], start = c(year(financeiras_naoexc[1,1]),
                                                               month(financeiras_naoexc[1,1])), freq = 12)
    }


    ## Período de construção do barômetro
    inicio = c(year(periodo[1]), month(periodo[1]))
    final = c(year(periodo[length(periodo)]), month(periodo[length(periodo)]))


    ## Base final no período de construção do barômetro
    covariaveisnovo_ts = window(covariaveis_ts, start = inicio, end = final)
    aux = apply(apply(covariaveisnovo_ts, 2, is.na), 2, sum)                  # Excluir séries com NA
    covariaveisnovo_ts = covariaveisnovo_ts[, names(aux[aux == 0])]           # Excluir séries com NA


    ## Teste de estacionariedade nas séries originais
    resultados_estac = data.frame(matrix(NA, ncol(covariaveisnovo_ts), 13))
    colnames(resultados_estac) = c("Serie", "KPSS_detrend", "KPSS_tau", "KPSS_mu",
                                   "ADF12_phi", "ADF24_phi", "ADF12_t", "ADF24_t",
                                   "DFGLS12_const", "DFGLS24_const", "DFGLS12_trend", "DFGLS24_trend",
                                   "Tend_det")

    for(j in 1:ncol(covariaveisnovo_ts)){

        serie = covariaveisnovo_ts[,colnames(covariaveisnovo_ts)[j]]
        resultados_estac[j,"Serie"] = colnames(covariaveisnovo_ts)[j]
        resultados_estac[j,"KPSS_detrend"] = teste_kpss(serie)$estac
        resultados_estac[j, "Tend_det"] = teste_kpss(serie)$tend_det
        resultados_estac[j,"KPSS_tau"] = teste_kpss2(serie, modelo = "tau")
        resultados_estac[j,"KPSS_mu"] = teste_kpss2(serie, modelo = "mu")
        resultados_estac[j,"ADF12_phi"] = teste_adf(serie, lag_max = 12)
        resultados_estac[j,"ADF24_phi"] = teste_adf(serie, lag_max = 24)
        resultados_estac[j,"ADF12_t"] = teste_adf2(serie, lag_max = 12)
        resultados_estac[j,"ADF24_t"] = teste_adf2(serie, lag_max = 12)
        resultados_estac[j,"DFGLS12_const"] = teste_dfgls(serie, modelo = "constant", lag_max = 12)
        resultados_estac[j,"DFGLS24_const"] = teste_dfgls(serie, modelo = "constant", lag_max = 24)
        resultados_estac[j,"DFGLS12_trend"] = teste_dfgls(serie, modelo = "trend", lag_max = 12)
        resultados_estac[j,"DFGLS24_trend"] = teste_dfgls(serie, modelo = "trend", lag_max = 24)

        # print(j)
    }


    ## Total de testes que concluíram pela NÃO estacionariedade da série
    resultados_estac[,"Soma_Resultados"] = rowSums(resultados_estac[,2:(ncol(resultados_estac) - 1)], na.rm = TRUE)

    # Caso o número de testes que concluíram pela NÃO estacionariedade de determinada série
    # seja igual ou superior a 5 concluiremos que a série NÃO é estacionária (Resultado Final = 1).
    # Caso contrário, diremos que a série É estacionária (Resultado Final = 0).
    resultados_estac[which(resultados_estac[,"Soma_Resultados"] >= 5),"Resultado Final"] = 1
    resultados_estac[which(resultados_estac[,"Soma_Resultados"] < 5),"Resultado Final"] = 0


    ## Séries estacionárias "em nível"
    series_E = resultados_estac[which(resultados_estac[,"Resultado Final"] == 0 &
                                          resultados_estac[,"Tend_det"] == 0 &
                                          !is.na(resultados_estac[,"Serie"])), "Serie"]
    if(length(series_E) != 0){
        estacionarias_originais = data.frame(Data = as.Date(covariaveisnovo_ts[,series_E]),
                                             covariaveisnovo_ts[,series_E])
    }else{
        estacionarias_originais = data.frame(Data = as.Date(covariaveisnovo_ts))
    }

    # Como o teste é feito nos resíduos do modelo pode ser que haja tendência determinística e o teste conclua
    # que a série é estacionária. Nesses casos devemos eliminar a tendência determinística.
    series_TD = resultados_estac[which(resultados_estac[,"Resultado Final"] == 0 &
                                           resultados_estac[,"Tend_det"] == 1 &
                                           !is.na(resultados_estac[,"Serie"])), "Serie"]

    if(length(series_TD) != 0){
        for(m in 1:length(series_TD)){

            # Obtendo os resíduos do modelo
            res = lm(covariaveisnovo_ts[,series_TD[m]] ~ seq(1, nrow(covariaveisnovo_ts), by = 1))$res
            estacionarias_originais = cbind(estacionarias_originais, res)
        }
        colnames(estacionarias_originais) = c("Data", series_E, series_TD)
    }


    ## Teste de estacionariedade das séries em primeira diferença (séries não estacionárias "em nível")
    series_NE = resultados_estac[which(resultados_estac[,"Resultado Final"] == 1 & !is.na(resultados_estac[,"Serie"])), "Serie"]

    if(length(series_NE) != 0){

        resultados_estac2 = data.frame(matrix(NA, length(series_NE), 13))
        colnames(resultados_estac2) = c("Serie", "KPSS_detrend", "KPSS_tau", "KPSS_mu",
                                        "ADF12_phi", "ADF24_phi", "ADF12_t", "ADF24_t",
                                        "DFGLS12_const", "DFGLS24_const", "DFGLS12_trend", "DFGLS24_trend",
                                        "Tend_det")
        for(k in 1:length(series_NE)){

            serie = covariaveisnovo_ts[,series_NE[k]]
            serie_dif = diff(serie, differences = 1)

            resultados_estac2[k,"Serie"] = series_NE[k]
            resultados_estac2[k,"KPSS_detrend"] = teste_kpss(serie_dif)$estac
            resultados_estac2[k,"Tend_det"] = teste_kpss(serie_dif)$tend_det
            resultados_estac2[k,"KPSS_tau"] = teste_kpss2(serie_dif, modelo = "tau")
            resultados_estac2[k,"KPSS_mu"] = teste_kpss2(serie_dif, modelo = "mu")
            resultados_estac2[k,"ADF12_phi"] = teste_adf(serie_dif, lag_max = 12)
            resultados_estac2[k,"ADF24_phi"] = teste_adf(serie_dif, lag_max = 24)
            resultados_estac2[k,"ADF12_t"] = teste_adf2(serie_dif, lag_max = 12)
            resultados_estac2[k,"ADF24_t"] = teste_adf2(serie_dif, lag_max = 12)
            resultados_estac2[k,"DFGLS12_const"] = teste_dfgls(serie_dif, modelo = "constant", lag_max = 12)
            resultados_estac2[k,"DFGLS24_const"] = teste_dfgls(serie_dif, modelo = "constant", lag_max = 24)
            resultados_estac2[k,"DFGLS12_trend"] = teste_dfgls(serie_dif, modelo = "trend", lag_max = 12)
            resultados_estac2[k,"DFGLS24_trend"] = teste_dfgls(serie_dif, modelo = "trend", lag_max = 24)

            # print(k)
        }


        ## Total de testes que concluíram pela NÃO estacionariedade da série
        resultados_estac2[,"Soma_Resultados"] = rowSums(resultados_estac2[,2:(ncol(resultados_estac2)-1)], na.rm = TRUE)

        # Caso o número de testes que concluíram pela NÃO estacionariedade de determinada série
        # seja igual ou superior a 5 concluiremos que a série NÃO é estacionária (Resultado Final = 1).
        # Caso contrário, diremos que a série É estacionária (Resultado Final = 0).
        resultados_estac2[which(resultados_estac2[,"Soma_Resultados"] >= 5),"Resultado Final"] = 1
        resultados_estac2[which(resultados_estac2[,"Soma_Resultados"] < 5),"Resultado Final"] = 0


        ## Séries estacionárias em primeira diferença
        seriesDIF_E = resultados_estac2[which(resultados_estac2[,"Resultado Final"] == 0 &
                                                  resultados_estac2[,"Tend_det"] == 0 &
                                                  !is.na(resultados_estac2[,"Serie"])), "Serie"]
        if(length(seriesDIF_E) != 0){
            estacionarias_DIF = data.frame(Data = as.Date(diff(covariaveisnovo_ts[,seriesDIF_E], differences = 1)),
                                           diff(covariaveisnovo_ts[,seriesDIF_E], differences = 1))
        }else{
            estacionarias_DIF = data.frame(Data = as.Date(diff(covariaveisnovo_ts)))
        }

        # Como o teste é feito nos resíduos do modelo pode ser que haja tendência determinística e o teste conclua
        # que a série é estacionária. Nesses casos devemos eliminar a tendência determinística.
        seriesDIF_TD = resultados_estac2[which(resultados_estac2[,"Resultado Final"] == 0 &
                                                   resultados_estac2[,"Tend_det"] == 1 &
                                                   !is.na(resultados_estac2[,"Serie"])), "Serie"]
        if(length(seriesDIF_TD) != 0){
            for(l in 1:length(seriesDIF_TD)){

                # Resíduos do modelo
                res = lm(covariaveisnovo_ts[-1,seriesDIF_TD[l]] ~ seq(1, (nrow(covariaveisnovo_ts)-1), by = 1))$res
                estacionarias_DIF = cbind(estacionarias_DIF, res)
            }
            colnames(estacionarias_DIF) = c("Data", seriesDIF_E, seriesDIF_TD)
        }
    }


    ## Matriz com as séries estacionárias (seja em nível ou em primeira diferença)
    covariaveis_final = merge(x = estacionarias_originais, y = estacionarias_DIF, by = "Data", all = TRUE)
    if(ncol(covariaveis_final) == 1){stop("Nenhuma variável é estacionária até a primeira diferença. Verifique as séries.")}


    ## Saídas da função
    if(!missing(sondagens)){

        sondagens_final = covariaveis_final[,which(substr(colnames(covariaveis_final), 1, 2) == "S_" |
                                                       colnames(covariaveis_final) == "Data")]

        if(!missing(financeiras)){

            financeiras_final = covariaveis_final[,which(substr(colnames(covariaveis_final), 1, 2) == "F_" |
                                                             colnames(covariaveis_final) == "Data")]
            saidas = list(covariaveis_estacionarias = covariaveis_final,
                          sondagens_estacionarias = sondagens_final,
                          financeiras_estacionarias = financeiras_final)
        }else{

            saidas = list(covariaveis_estacionarias = covariaveis_final,
                          sondagens_estacionarias = sondagens_final)
        }
    }else{

        financeiras_final = covariaveis_final[,which(substr(colnames(covariaveis_final), 1, 2) == "F_" |
                                                         colnames(covariaveis_final) == "Data")]
        saidas = list(covariaveis_estacionarias = covariaveis_final,
                      financeiras_estacionarias = financeiras_final)
    }

    saidas = append(saidas,
                    list(resultados = resultados_estac, resultados_dif = resultados_estac2))
    return(saidas)
}
