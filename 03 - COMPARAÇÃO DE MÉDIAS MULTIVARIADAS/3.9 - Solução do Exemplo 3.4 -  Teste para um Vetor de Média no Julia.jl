#Teste para um vetor de médias
#SOLUÇÃO DO EXEMPLO 3.4 NO JULIA
using DataFrames, Statistics, LinearAlgebra
#Matriz de dados para o volume sem casca e volume com casca
X=[7.63 9.50
7.90 9.80
7.98 9.88
8.07 9.98
8.17 10.13
8.32 9.45
8.33 10.30
8.35 10.31
8.41 10.26
8.53 10.13]
#Vetor de médias a ser comparado (valores esperados)
mizero=[8 10]
#Vetor de médias
Media=mean(X, dims=1)
#Matriz variância-covariância
S=cov(X)
75
#Inversa da matriz variância-covariância
invS=inv(S)
#Tamanho da amostra (n) e número de variáveis (p)
(n,p)=size(X)
#Cálculo da estatística T2 de Hotelling
T2=n*(Media-mizero)*invS*(Media-mizero)'
#Cálculo do valor crítico de F
F = (((n - 1) * p) / (n - p)) * quantile(FDist(p, n - p), 0.95)