#Teste para vetores independentes com matrizes variância-covariância heterogêneas
#SOLUÇÃO DO EXEMPLO 3.7 NO JULIA
using DataFrames, Statistics, LinearAlgebra
#Matriz de dados para o grupo 1 (Preparo do solo A)
X1=[0.6 0.3 33 11 0.09
0.6 0.3 64 7 0.04
0.2 0.3 71 6 0.06
0.1 0.7 50 9 0.05
0.6 0.4 98 6 0.07
0.2 0.6 24 9 0.08
0.2 0.4 66 5 0.06
0.4 0.4 105 8 0.08
0.8 0.6 113.8 7 0.1
0.6 0.4 73 6 0.05
0.9 0.4 78 5 0.07
0.1 0.5 31 10 0.05
0.6 0.5 171 8 0.08
0.2 0.6 78 9 0.06
0.2 0.4 56 6 0.06
0.4 0.4 78 8 0.07
0.1 0.3 64 6 0.04
0.2 0.4 50 10 0.07
0.6 0.4 24 8 0.07
0.1 0.4 32 5 0.07
0.2 0.5 71 6 0.08
0.4 0.4 76 6 0.06
0.2 0.4 18 7 0.08
0.2 0.6 29 6 0.07
0.4 0.4 29 8 0.1
0.4 0.3 73 6 0.07
0.4 0.4 71 13 0.09
0.7 0.5 146 7 0.08
0.2 0.5 20 7 0.12
0.2 0.6 91 11 0.05
0.1 0.3 93 8 0.1
0.2 0.6 21 10 0.12
0.6 0.5 91 6 0.08
0.2 0.3 53 7 0.07
0.4 0.4 54 7 0.06
0.4 0.4 49 4 0.1
0.2 0.4 71 8 0.07]
#Matriz de dados para o grupo 2 (Preparo do solo B)
X2=[0.6 0.3 127 10 0.08
0.6 0.3 87 7 0.06
0.6 0.3 69 20 0.07
85
0.2 0.4 58 8 0.05
0.6 0.3 42 8 0.04
0.6 0.4 32 7 0.04
0.2 0.3 37 8 0.08
0.8 0.5 52 20 0.07
0.9 0.4 87 5 0.07
0.6 0.5 44 8 0.06
0.4 0.4 49 8 0.06
0.4 0.7 61 15 0.12
0.2 0.3 78 7 0.07
0.2 0.4 53 11 0.07
1.1 0.5 75 8 0.09
0.9 0.4 73 5 0.07
0.6 0.5 99 6 0.07
0.4 0.4 78 9 0.06
0.4 0.3 40 6 0.1
0.9 0.3 80 5 0.06
0.8 0.3 76 13 0.06
0.9 0.3 87 10 0.08
0.6 0.3 44 10 0.07
0.8 0.3 171 10 0.06
0.9 0.5 40 9 0.08
0.6 0.4 57 9 0.09
0.9 0.4 104 7 0.07
0.8 0.4 71 7 0.05
0.8 0.3 66 6 0.05
1.2 0.4 118 6 0.08
0.4 0.4 41 8 0.13
0.9 0.5 107 12 0.6
0.1 0.4 41 7 0.1
0.4 0.6 95 20 0.12
0.2 0.8 79 13 0.07
0.4 0.5 59 8 0.09
0.4 0.4 82 6 0.05]
#Tamanho da amostra do grupo 1 (n1) e número de variáveis (p)
(n1,p)=size(X1) 
#Tamanho da amostra do grupo 2 (n2) e número de variáveis (p)
(n2,p)=size(X2) 
#Vetor de médias para o grupo 1
MEDIAX1=mean(X1,dims=1)
#Vetor de médias para o grupo 2
MEDIAX2=mean(X2,dims=1)
#Matriz variância-covariância para o grupo 1
S1=cov(X1)
86
#Matriz variância-covariância para o grupo 2
S2=cov(X2)
#Matriz variância-covariância ponderada
SP=((1/n1)*S1)+((1/n2)*S2)
#Cálculo da estatística T2 de Hotelling
T2=(MEDIAX1-MEDIAX2)*inv(SP)*(MEDIAX1-MEDIAX2)'
#Cálculo do valor de Qui-quadrado crítico
QUI2 = quantile(Chisq(p), 0.95)