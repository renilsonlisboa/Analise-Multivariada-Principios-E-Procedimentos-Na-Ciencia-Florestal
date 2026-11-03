#Análise da homogeneidade de matrizes variância-covariância
#SOLUÇÃO DO EXEMPLO 3.3 NO JULIA
using DataFrames, Statistics, LinearAlgebra
#Matriz de dados do grupo 1 (clone A)
XA=[0.05627 0.04566
0.18012 0.15376
0.29275 0.24932
0.05525 0.04569
0.12055 0.10179
0.21626 0.18614
0.05468 0.04580
0.11396 0.09670
0.25274 0.20748
0.04135 0.03316
0.11727 0.09503
0.14707 0.11882
0.03518 0.03071
0.13309 0.11207
0.22043 0.18568
0.04550 0.03724
0.09090 0.07592
0.15289 0.12450
0.04923 0.03961
0.13854 0.11848
0.22338 0.19068
0.04144 0.03249
0.21322 0.18033
0.18346 0.15311
0.04857 0.04109
0.13621 0.11484
0.31096 0.26852
0.05460 0.04813
0.12250 0.09829
0.19152 0.16087]
#Matriz de dados do grupo 2 (clone B)
XB=[0.04058 0.03349
0.09321 0.07654
0.17324 0.14701
0.15231 0.13220
72
0.20407 0.17811
0.07278 0.06258
0.06487 0.05284
0.19848 0.16709
0.15004 0.12780
0.03816 0.03125
0.11024 0.09459
0.18531 0.15747
0.04311 0.03596
0.15980 0.12702
0.33257 0.28749
0.18213 0.15175
0.38640 0.33046
0.16669 0.13862
0.01863 0.01478
0.06900 0.05616
0.12995 0.10654
0.05183 0.04243
0.13582 0.11089
0.23464 0.19872
0.03010 0.02526
0.11272 0.09614
0.21245 0.17713
0.04858 0.03935
0.11949 0.09769
0.23579 0.20158]
#Matriz variância-covariância do grupo1
S1=cov(XA)
#Matriz variância-covariância do grupo 2
S2=cov(XB)
#Número de grupos
k=2
#Tamanho da amostra do grupo 1 (n1) e número de variáveis (p)
(n1,p)=size(XA)
#Tamanho da amostra do grupo 2 (n2) e número de variáveis (p)
(n2,p)=size(XB)
#Tamanho total das amostras
n=n1+n2
#Graus de liberdade para o grupo 1
glX1=n1-1
#Graus de liberdade para o grupo 2
glX2=n2-1
73
#Matriz variância-covariância combinada
S=(1/(glX1+glX2))*((glX1*S1)+(glX2*S2))
#Determinante da matriz variância-covariância combinada
DETS=det(S)
#Determinante da matriz variância-covariância para o grupo 1
DETS1=det(S1)
#Determinante da matriz variância-covariância combinada para o grupo 2
DETS2=det(S2)
#Cálculo da estatística M de Box
M=(glX1*log(DETS)+glX2*log(DETS))-(glX1*log(DETS1)+glX2*log(DETS2))
#Fator de escala para amostras de mesmo tamanho
C1=1-(2*p^2+3*p-1*k+1)/(6*(p+1)*k*n)
#Cálculo de M de Box corrigido
MC1=M*C1
#Graus de liberdade para obtenção do valor de Qui-Quadrado
Gl=0.5*(k-1)*p*(p+1)
#Valor crítico de Qui-Quadrado 
Qui2 = quantile(Chisq(Gl), 0.95)