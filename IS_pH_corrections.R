devtools::document()
devtools::load_all()

####EDTA####

IS = 0.5*(0.24+0.24+0.14+0.14+0.01+0.01)

?log10K.0.calculator

log10K.0 = log10K.0.calculator(8.79, 0.1, 2, -4, -2)

log10K.0

?log10K.IS.calculator

log10.K.IS =  log10K.IS.calculator(8.79, 0.1, 0.1, 2, -4, -2)

log10.K.IS

log10.K.IS =  log10K.IS.calculator(8.79, 0.1, IS, 2, -4, -2)

log10.K.IS

x = c(1:1500)
y = c()

for (i in x){
  y[i] = log10K.IS.calculator(8.79, 0.1, x[i]/1000, 2, -4, -2)
}

plot(x/1000, y)

?pKa.IS.calculator

pKa1 = pKa.IS.calculator(10.17, IS, 0.1, -4, -3)

pKa2 =  pKa.IS.calculator(6.11, IS, 0.1, -3, -2)

a1 = 1/(1 + 10^(pKa1 - 7) + 10^(pKa2 + pKa1 - 14))

a1

K.app = a1*10^log10.K.IS

1000/K.app

####ATP####

IS = 0.5*(0.24+0.24+0.14+0.14+0.01+0.01)

?log10K.0.calculator

log10K.0 = log10K.0.calculator(4.3, 0.1, 2, -4, -2)

log10K.0

?log10K.IS.calculator

log10.K.IS =  log10K.IS.calculator(4.3, 0.1, 0.1, 2, -4, -2)

log10.K.IS

log10.K.IS =  log10K.IS.calculator(4.3, 0.1, IS, 2, -4, -2)

log10.K.IS

x = c(1:1500)
y = c()

for (i in x){
  y[i] = log10K.IS.calculator(4.3, 0.1, x[i]/1000, 2, -4, -2)
}

plot(x, y)

?pKa.IS.calculator

pKa1 = pKa.IS.calculator(6.53, IS, 0.1, -4, -3)

pKa2 =  pKa.IS.calculator(4, IS, 0.1, -3, -2)

a1 = 1/(1 + 10^(pKa1 - 7) + 10^(pKa2 + pKa1 - 14))

a1

K.app = a1*10^log10.K.IS

1000/K.app

####Glu 6P####

####Glutamate####

IS = 0.5*(0.24+0.24+0.14+0.14+0.01+0.01+0.1+0.1)

?log10K.0.calculator

log10K.0 = log10K.0.calculator(1.9, 0.1, 2, -2, 0)

log10K.0

?log10K.IS.calculator

log10.K.IS =  log10K.IS.calculator(1.9, 0.1, 0.1, 2, -2, 0)

log10.K.IS

log10.K.IS =  log10K.IS.calculator(1.9, 0.1, IS, 2, -2, -0)

log10.K.IS

x = c(1:1500)
y = c()

for (i in x){
  y[i] = log10K.IS.calculator(4.3, 0.1, x[i]/1000, 2, -4, -2)
}

plot(x, y)

?pKa.IS.calculator

pKa1 = pKa.IS.calculator(9.59, IS, 0.1, -2, -1)

pKa2 =  pKa.IS.calculator(4, IS, 0.1, -3, -2)

a1 = 1/(1 + 10^(pKa1 - 7) + 10^(pKa2 + pKa1 - 14))

a1

K.app = a1*10^log10.K.IS

1000/K.app
