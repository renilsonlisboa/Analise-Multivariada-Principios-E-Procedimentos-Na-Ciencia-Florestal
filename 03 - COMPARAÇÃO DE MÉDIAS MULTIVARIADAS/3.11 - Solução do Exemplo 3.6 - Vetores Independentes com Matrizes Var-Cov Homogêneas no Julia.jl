#Teste para vetores independentes com matrizes variância-covariância homogêneas
#SOLUÇÃO DO EXEMPLO 3.6 NO JULIA
using DataFrames, Statistics, LinearAlgebra
#Matriz de dados para o grupo 1 (Componente copa)
X1=[72.25 25.02 2.73
71.54 25.67 2.78
71.79 25.56 2.65
73.21 25.06 1.73
70.64 26.74 2.62
70.70 27.12 2.19
69.65 27.39 2.96
69.86 28.04 2.10
75.71 23.50 0.79
72.84 25.09 2.07
71.20 26.13 2.67
70.48 27.44 2.08
71.76 25.36 2.88
69.26 28.09 2.65
72.51 25.89 1.61
69.42 28.84 1.74
71.80 25.89 2.30]
#Matriz de dados para o grupo 2 (Componente Fuste)
X2=[77.71 21.49 0.79
77.04 22.02 0.94
79.34 19.95 0.72
79.02 19.26 1.72
77.05 22.35 0.60
75.09 23.82 1.09
77.86 21.19 0.95
76.30 22.64 1.06
79.82 18.19 1.99
74.61 23.53 1.86
79.59 19.70 0.70
78.31 20.79 0.89
75.24 23.39 1.36
75.26 23.19 1.55
78.08 21.14 0.78
76.54 21.65 1.81
79.34 20.03 0.63]
#Tamanho da amostra do grupo 1 (n1) e número de variáveis (p)
(n1,p)=size(X1) 
#Tamanho da amostra do grupo 2 (n2) e número de variáveis (p)
(n2,p)=size(X2) 
81
#Vetor de médias para o grupo 1
MEDIAX1=mean(X1,dims=1)
#Vetor de médias para o grupo 2
MEDIAX2=mean(X2,dims=1)
#Matriz variância-covariância para o grupo 1
S1=cov(X1)
#Matriz variância-covariância para o grupo 2
S2=cov(X2)
#Matriz variância-covariância combinada
Sp=((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
#Cálculo da estatística T2 de Hotelling
T2=(MEDIAX1-MEDIAX2)*(inv((1/n1 + 1/n2)*(Sp)))*(MEDIAX1-MEDIAX2)'
#Transformação de T2 de Hotelling para o valor de F calculado
Fc=(T2*(n1+n2-p-1))/((n1+n2-2)*p)
#Cálculo do valor de F crítico
FT = quantile(FDist(p, n1 + n2 - p - 1), 0.95)