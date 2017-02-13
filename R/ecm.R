#' An Error Correction Model Function
#'
#' This function transforms your OLS coefficients into error correction coefficients. It also returns a list with the p-values for stargazer/texreg packages.
#' @param formula A regression formula.
#' @param serie A dataframe or ts object with the series used by the formula provided.
#' @keywords Error Correction Model, ECM
#' @export
#' @import zoo dynlm
#' @author Fernando Teixeira
#' @examples
#' ECM(y ~ x + z, df)

ECM = function (formula, serie) {
    
    assign("formula1", formula, envir = .GlobalEnv)
    
    for (j in 1:length(attributes(serie)$class)) {
        if (attributes(serie)$class[j] != "mts" | attributes(serie)$class[j] != 
            "ts") {
            tserie = ts(serie)
            serie = tserie
        }
    }
    assign("tserie", serie, envir = .GlobalEnv)
    
    a = dynlm(formula1, data = tserie)
    f = as.character(formula1)
    f = f[2]
    f = strsplit(f, split = ",")
    f0 = substr(f[[1]][1], start = 1, stop = 2)
    if (f0 == "d(") {
        f1 = substr(f[[1]][1], start = 3, stop = 100)
    } else {
        f1 = substr(f[[1]][1], start = 6, stop = 100)
    }
    nomes = names(a$coefficients)
    nome = NULL
    for (i in 1:length(nomes)) {
        nome[i] = strsplit(nomes[i], split = ",")
        f3 = substr(nome[[i]][1], start = 3, stop = 100)
        if (f1 == f3) {
            f4 = nomes[i]
            indice = i
        }
    }
    old.coef = a$coefficients
    su = summary(a)
    pval = list(su$coefficients[, 4])
    for (j in 1:length(nomes)) {
        if (j != indice) {
            f2 = substr(nome[[j]][1], start = 1, stop = 2)
            if (f2 == "L(") {
                su$coefficients[j, 1] = -su$coefficients[j, 1]/a$coefficients[f4]
                a$coefficients[j] = -a$coefficients[j]/a$coefficients[f4]
            }
        }
    }
    lista = list(new.summary = su, new.coef = a, pval = pval, 
                 old.coef = old.coef)
    return(invisible(lista))
}