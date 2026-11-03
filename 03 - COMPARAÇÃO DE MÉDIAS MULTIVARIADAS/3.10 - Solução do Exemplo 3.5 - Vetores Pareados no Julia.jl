#Teste para vetores pareados
#SOLUÇÃO DO EXEMPLO 3.5 NO JULIA
using DataFrames, Statistics, LinearAlgebra
#Matriz de dados para o grupo 1 (Método A)
X1=[50.2 2.9 2.1
54.1 3.1 2.2
48.9 2.4 2.4
50.6 2.7 2.3
57.4 3.0 1.9
49.6 3.7 1.8
45.7 2.8 2.0
77
44.0 2.9 2.0
52.1 3.2 2.2
46.7 3.1 1.7
49.3 3.3 1.9
51.8 2.8 2.0]
#Matriz de dados para o grupo 2 (Método B)
X2=[49.9 2.8 1.9
52.3 3.2 2.3
46.7 2.2 2.5
49.0 2.9 2.0
54.1 3.0 1.7
48.2 3.2 2.0
46.4 3.0 1.8
44.0 2.7 1.6
49.6 2.9 2.0
44.7 3.3 1.5
50.1 3.4 1.9
49.9 2.9 2.4]
#Matriz das diferenças entre grupos
DIF=X1-X2
#Matriz variância-covariância das diferenças
SDIF=cov(DIF)
#Inversa da matriz variância-covariância das diferenças
INVSDIF=inv(SDIF)
#Médias das diferenças
DIFMEDIA=mean(DIF,dims=1)
#Tamanho da amostra (n) e número de variáveis (p)
(n,p)=size(DIF)
#Cálculo da estatística T2 de Hotelling
T2=n*DIFMEDIA*INVSDIF*DIFMEDIA'
#Cálculo do valor de F crítico corrigido
F = ((n - 1) * p) / (n - p) * quantile(FDist(p, n - p), 0.99)
