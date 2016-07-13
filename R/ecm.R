#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

ECM = function(formula, serie){
  
  
  for (j in 1: length(attributes(serie)$class)){
    if (attributes(serie)$class[j] != "mts" | attributes(serie)$class[j] != "ts"){
      tserie = ts(serie)
      serie = tserie
    }   
    
  }
  
  assign("tserie", serie, envir = .GlobalEnv)
  
  formula = formula 
  a=dynlm(formula, data = tserie)
  f = as.character(formula)
  f = f[2]
  f = strsplit(f, split = ",")
  f1 = substr(f[[1]][1], start = 6, stop = 100)
  nomes = names(a$coefficients)
  nome = NULL
  
  for (i in 1:length(nomes)){
    
    nome[i] = strsplit(nomes[i], split = ",")
    f3 = substr(nome[[i]][1], start=3, stop = 100)    
    if (f1 == f3){
        f4 = nomes[i]
        index = i
    }
  }

  for (j in 1: length(nomes)){
      if (j != index){
          f2 = substr(nome[[j]][1], start=1, stop = 2)
          if(f2 == "L("){
            a$coefficients[j] = - a$coefficients[j]/a$coefficients[f4]              
          }
      }
  }
  
  
  lista = list('ecm' = summary(a), 'coefficients' = a)
  return(invisible(lista))
}