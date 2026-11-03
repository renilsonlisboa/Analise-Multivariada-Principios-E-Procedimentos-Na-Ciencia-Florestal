# SOLUÇÃO DO EXEMPLO 6.6 NO JULIA
#Matriz das variáveis respostas
Y=[
 43.48 8.42 73.70 0.56 0.49
 67.96 12.79 121.04 0.56 0.50
 65.19 11.61 114.27 0.54 0.47
185
 70.50 14.36 119.43 0.59 0.50
 47.40 8.53 78.21 0.59 0.51
 52.90 10.58 95.22 0.57 0.49
 61.30 12.87 109.73 0.53 0.46
 59.70 11.04 98.51 0.57 0.48
 72.10 17.46 143.16 0.52 0.53
 70.87 20.93 138.80 0.51 0.50
 85.31 22.18 197.15 0.51 0.50
 66.04 20.70 124.67 0.52 0.54
 74.06 17.77 216.39 0.51 0.51
 89.53 18.80 222.53 0.54 0.54
 69.15 18.58 119.72 0.58 0.55
 87.99 21.13 146.65 0.58 0.56
 110.08 19.43 178.07 0.58 0.58
 159.37 33.56 278.88 0.58 0.57
 127.09 17.30 222.87 0.60 0.56
 111.52 14.38 201.16 0.58 0.54
 143.13 24.24 257.03 0.56 0.56
 137.43 27.21 236.47 0.58 0.58
 115.42 20.30 190.76 0.58 0.57
 152.91 24.84 246.90 0.58 0.53]
#Matriz de variáveis categóricas para representar os grupos (tratamentos) 
X=[ 
 1 1 0 0
 1 1 0 0
 1 1 0 0
 1 1 0 0
 1 1 0 0
 1 1 0 0
 1 1 0 0
 1 1 0 0
 1 0 1 0
 1 0 1 0
 1 0 1 0
 1 0 1 0
 1 0 1 0
 1 0 1 0
 1 0 1 0
 1 0 1 0
 1 0 0 1
 1 0 0 1
 1 0 0 1
 1 0 0 1
 1 0 0 1
 1 0 0 1
 1 0 0 1
 1 0 0 1]
#Determinação do tamanho da matriz X
(m,n)=size(X)
186
#Somas de quadrados e produtos não corrigidas
S=X'X
#Matriz dos coeficientes / inversa de X'X obtida pelo método de moore penrose
B=pinv(X'X)X'Y
#Matriz de somas de quadrados dos resíduos corrigida
SQRes=Y'Y-B'X'Y
#Matriz de somas de quadrados dos tratamentos não corrigida
SQTrat_NC=B'X'Y
#Matriz de soma de quadrados totais não corrigida
SQTot_NC=Y'Y
#Somas de quadrados corrigidas
#Matriz de desvios
a=Y-((ones(m)*ones(m)')*Y).*(1/m)
#Matriz de somas de quadrados totais corrigida
SQTot_Corr=a'*a
#Matriz de somas de quadrados dos tratamentos corrigida
SQTrat_Corr=SQTot_Corr-SQRes
#Cálculo de Lambda de Wilks
L=det(SQRes)/det(SQTrat_Corr+SQRes)
#Número de variáveis
p=size(Y,2)
#Número de grupos
g=n-1
#Cálculo do valor de Lambda de Wilks corrigido
Lcorrigido=-(m-1-((p+g)/2))*log(L)
#Valor crítico de Qui-Quadrado 
QUIQUADRAINV = quantile(Chisq(p * (g - 1)), 0.99)