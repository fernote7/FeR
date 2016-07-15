#' ADF Function
#'
#' This function allows you to run an ADF test through several contiguous columns of a dataframe using urca package.
#' @param data Dataframe.
#' @param col.init Initial column.
#' @param col.fin Final column.
#' @param lags maximum lags for test. Defaults to 24.
#' @param select.lags Information criterion ("Fixed", "AIC", "BIC"). Defaults to AIC.
#' @param mymfrow Plot parameters. Defaults to 5 lines and 3 colunns of plots.
#' @param mymar Plot parameters. Size of plot margins, defaults to c(3,2,1.5,1).
#' @param myoma Plot parameters. Size of plot outside margins, defaults to c(0, 0, 2, 0).
#' @param dif Apply test to first difference? Defaults to 'y'.
#' @param res Compute residuals? Defaults to 'y'.
#' @param precision Number of digits in a number.
#' @param plt Plot? Defaults to FALSE.
#' @keywords ADF test urca
#' @export
#' @examples
#' func_adf(series,1,3, dif = 'y', res = 'n', mymfrow = c(3,1), lags=13)
#' @author Fernando Teixeira

func_adf <- function(data, col.init, col.fin, lags = 24, select.lags = "AIC", 
                     mymfrow = c(5,3), mymar = c(3,2,1.5,1), myoma = c(0, 0, 2, 0), 
                     dif='y', res='y', precision = 2, plt=F){
    
    if(dif != 'y' && dif!= 'n' || res!= 'y' && res!= 'n'){
        
        stop("Diferenciar assume valor y ou n")
    }
    
    
    numcol = col.fin-col.init + 1 # núm total de colunas
    b = rep(1, numcol)
    c = rep(1,3)
    c1 = rep(1, numcol)
    d = matrix(nrow = length(data[,1])-lags-1, ncol= numcol)
    residuos = matrix(nrow = length(data[,1])-lags-1, ncol= numcol)
    residuos2 = matrix(nrow = length(data[,1])-lags-2, ncol= numcol)
    f = rep(NA,numcol)
    trend = rep(NA,numcol)
    ttrend = rep(NA,numcol)
    select.lags = select.lags
    nomescol=colnames(data[,col.init:col.fin])
    par(mfrow=mymfrow, mar=mymar, oma=myoma)
    precision = precision
    
    
    for (i in 1:numcol){
        a = ur.df(data[,col.init+i-1], lags = lags, 
                  selectlags = select.lags, type=c("trend"))
        sumario=summary(a)
        residuos[,i] = a@res
        trend[i] ="trend e constante"
        if (sumario@testreg$coefficients[3,4] > 0.05){
            trend[i] ="constante"
            a = ur.df(data[,col.init+i-1], lags = lags, 
                      selectlags = select.lags, type=c("drift"))
            sumario=summary(a)
            residuos[,i] = a@res
        
            if (sumario@testreg$coefficients[1,4] > 0.05){
            trend[i] ="sem trend ou constante"
            a = ur.df(data[,col.init+i-1], lags = lags, 
                      selectlags = select.lags, type=c("none"))
            residuos[,i] = a@res
            }
        }
        b[i] = a@teststat[1]
        c = a@cval
        c1[i] = a@cval[1,2]
        d[,i] = a@res
        
        if (plt == TRUE){
            acf(a@res, lag.max = 36, drop.lag.0= T, main=paste0("ACF",i))
        }    
        
        if (a@teststat[1]<c[1,2]){
            f[i] = 'série estacionária'
        }
        else{
            f[i] = 'série não-estacionária'
        }
    }

    if (plt == TRUE){    
        mtext('Funções de Autocorrelação', outer = TRUE, cex = 1.5)
    }
    
    par(mfrow=mymfrow, mar=mymar, oma=myoma)
    
    if (plt == TRUE){
        if(res=='y'){
            for (i in 1:numcol){
                
                matplot(residuos[,i], type = "l", pch = 20, main=paste0("Resíduos",i))
                
            }
            mtext('Gráfico Resíduos', outer = TRUE, cex = 1.5)
        }
    }
        
    if (dif == 'y'){
        par(mfrow=mymfrow, mar=mymar, oma=myoma)
        bb= rep(1, numcol)
        d2 = matrix(nrow = length(data[,1])-lags-2, ncol= numcol)
        f2 = rep(NA,length(grep('série não-estacionária', f)))
        nomescol2 = rep(NA,numcol)
        cc = rep(1,3)
        c2 = rep(1, numcol)
        serie_diff = matrix(nrow = length(data[,1])-1, ncol= numcol)
        
        for (i in 1:length(f)){ 
            if(f[i] == 'série não-estacionária'){
                aa = ur.df(diff(data[,col.init+i-1], k=1, differences = 1), 
                           lags = lags, selectlags = select.lags, type=c("trend"))
                sumario=summary(aa)
                residuos2[,i] = aa@res
                serie_diff[,i] = diff(data[,col.init+i-1], k=1, differences = 1)
                ttrend[i] ="trend e constante"
                if (sumario@testreg$coefficients[3,4] > 0.05){
                    ttrend[i] ="constante"
                    aa = ur.df(diff(data[,col.init+i-1], k=1, differences = 1), 
                               lags= lags, selectlags = select.lags, type=c("drift"))
                    sumario=summary(aa)
                    residuos2[,i] = aa@res
                
                    if (sumario@testreg$coefficients[1,4] > 0.05){
                    ttrend[i] ="sem trend ou constante"
                    aa = ur.df(diff(data[,col.init+i-1], k=1, differences = 1), 
                               lags = lags, 
                               selectlags = select.lags, type=c("none"))
                    residuos2[,i] = aa@res
                    }
                }
                bb[i] = format(round(aa@teststat[1] ,precision), nsmall = precision)
                cc = aa@cval
                
                d2[,i] = aa@res
                if(f[i] == 'série não-estacionária'){
                    if (plt == TRUE){
                        acf(aa@res, lag.max = 36, drop.lag.0= T, main=paste0("ACF diff",i))
                    }
                    c2[i] = aa@cval[1,2]
                    bb[i] = format(round(aa@teststat[1] ,precision), 
                                   nsmall = precision)
                    
                }
                if (aa@teststat[1]<cc[1,2]){
                    f2[i] = 'série estacionária'
                    c2[i] = aa@cval[1,2]
                    bb[i]=format(round(aa@teststat[1] ,precision), nsmall = precision)
                    
                }
                else{
                    f2[i] = 'série não-estacionária'
                    c2[i] = aa@cval[1,2]
                    bb[i]=format(round(a@teststat[1], precision), nsmall = precision)
                }
            } 
            else{  
                f2[i] = 'série já era estacionária'
                ttrend[i] = "-"
                c2[i] = "-"
                bb[i] = "-"
            }
            
        }
        if (plt == TRUE){
            mtext('Funções de Autocorrelação diff', outer = TRUE, cex = 1.5)
        }
        
        par(mfrow=mymfrow, mar=mymar, oma=myoma)
        
        if (plt == TRUE){
            if(res == 'y'){
                for (i in 1:numcol){
                    if (!is.na(residuos2[1,i])){  
                        matplot(residuos2[,i], type = "l", pch = 20, 
                                main=paste0("Resíduos diff",i))
                    }
                }
                mtext('Gráfico Resíduos das Séries diff.', outer = TRUE, cex = 1.5)
            }
        }
        
        if (plt == TRUE){
            par(mfrow=mymfrow, mar=mymar, oma=myoma)
            for (i in 1:numcol){
                if (!is.na(residuos2[1,i])){
                    matplot(serie_diff[,i], type = "l", 
                            main=paste0("Série diff",i))
                }
            }
        
        mtext('Gráfico Séries Diferenciadas', outer = TRUE, cex = 1.5)
        }
    }
    
    suppressWarnings(par())
    
    if(dif == 'y'){
        nomescol=colnames(data[,col.init:col.fin])
        b = as.numeric(b)
        b = format(round(b, precision), nsmall = precision)
        nomescol = as.data.frame(cbind(f,b,c1,trend,f2,bb, c2, ttrend), 
                                 row.names = nomescol)
        colnames(nomescol) = c("Teste ADF", "Estatística de teste", 
                               "Valor Crítico","Tipo teste 1",
                               "Teste ADF diff", "Estatística de teste 2", 
                               "Valor Crítico", "Tipo teste 2")
        lista=list('nomescol'=nomescol,'serie_diff'=serie_diff)
    }
    else if (dif == 'n') {
        nomescol=colnames(data[,col.init:col.fin])
        b=format(round(b, precision), nsmall = precision)
        nomescol = as.data.frame(cbind(f,b,c1,trend), row.names = nomescol)
        colnames(nomescol) = c("Teste ADF", "Estatística de teste", 
                               "Valor Crítico","Tipo teste 1")
        lista=list('nomescol'= nomescol)
    }
    
    return(lista)
}