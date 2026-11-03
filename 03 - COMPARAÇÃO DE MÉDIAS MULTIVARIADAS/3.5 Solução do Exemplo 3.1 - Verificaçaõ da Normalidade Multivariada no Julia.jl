#SOLUÇÃO DO EXEMPLO 3.1 NO JULIA
using DataFrames, CSV, Statistics, LinearAlgebra, Distributions, PyPlot,
LinearRegression
#Variáveis DAP e Altura
X=[19.58 19.7
7.48 14.9
16.07 18.2
9.23 15.7
9.39 11.8
5.41 12
8.44 12
10.19 18.5
18.94 18.8
23.55 21.7]
#Média
MEDIA=mean(X, dims=1)
#Matriz variância-covariância
S=cov(X)
#Inversa da matriz variância-covariância
INVS=inv(S)
62
#Cálculo da distância generalizada
D2j=sum(((X.-MEDIA)*INVS).*(X.-MEDIA),dims=2)
#Tamanho da amostra (n) e número de variáveis (p)
(n,p)=size(X)
#Obtenção do valor de Qui-quadrado tabelado
QUIQUADTAB = quantile(Chisq(p), 0.5)
#Filtro de valores menores ou igual ao valor de Qui-Quadrado tabelado
Contagem=(filtered = D2j[(D2j[:,:] .<= VALORDEQUI)])
#Número de valores menores ou iguais ao valor de Qui-Quadrado tabelado
ProporçãoMenor=length(Contagem)