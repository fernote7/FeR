serie = cbind(rnorm(100),rnorm(100))
serie = as.data.frame(serie)
formula = diff(V1) ~ diff(V2) + L(V1) + L(V2)
serie = ts(serie)

t1 = dynlm(formula, serie)
t = summary(t1)

# Extract the p-values
pvals <- coef(t)[,4]

# Use the symnum function to produce the symbols
sigSymbols <- symnum(pvals, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("$^{***}$", "$^{**}$", "$^{*}$", "$.$", ""))



t$coefficients[,4] = paste0(t$coefficients[,4], sigSymbols)
t2=xtable(t)


print(xtable(t, digits = -3), type='latex', sanitize.text.function=identity, digits=2) 
