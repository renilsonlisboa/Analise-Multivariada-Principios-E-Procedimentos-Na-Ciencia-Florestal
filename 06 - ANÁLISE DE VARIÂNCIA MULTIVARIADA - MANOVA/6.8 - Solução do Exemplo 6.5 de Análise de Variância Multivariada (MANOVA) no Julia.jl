#MÉTODO 1: OBTENÇÃO DA MANOVA POR MEIO DO CÁLCULO DOS DESVIOS
using DataFrames, CSV, Statistics, LinearAlgebra
#Matriz de dados geral
X=[0.50 69.33 27.43 3.24 5113.53
0.53 69.01 27.81 3.18 5115.39
0.51 68.94 27.88 3.19 5043.65
0.50 69.07 27.61 3.32 5081.57
0.51 75.24 22.58 2.18 4708.60
0.51 73.24 24.28 2.49 4553.28
0.51 73.89 23.89 2.21 4479.64
0.48 73.93 23.75 2.33 4664.67
0.52 73.44 23.91 2.65 4754.10
0.54 72.15 24.97 2.88 4909.00
0.50 73.81 23.48 2.71 4704.58
0.49 73.80 23.71 2.49 4690.54
0.55 74.63 23.50 1.88 4620.83
0.61 73.63 24.13 2.24 4647.88
0.60 72.24 25.23 2.53 4657.54
0.59 72.69 24.87 2.45 4693.50]
#Matriz de dados para o grupo 1 (Tratamento 1)
X1=X[1:4,:]
#Matriz de dados para o grupo 2 (Tratamento 2)
X2=X[5:8,:]
177
#Matriz de dados para o grupo 3 (Tratamento 3)
X3=X[9:12,:]
#Matriz de dados para o grupo 4 (Tratamento 4)
X4=X[13:16,:]
#Tamanho da amostra do grupo 1 (n1) e número de variáveis (p)
(n1,p)=size(X1)
#Tamanho da amostra do grupo 2 (n2) e número de variáveis (p)
(n2,p)=size(X2)
#Tamanho da amostra do grupo 3 (n3) e número de variáveis (p)
(n3,p)=size(X3)
#Tamanho da amostra do grupo 4 (n4) e número de variáveis (p)
(n4,p)=size(X4)
#Tamanho total das amostras
n=n1+n2+n3+n4
#Vetor do tamanho das amostras dos grupos
ng=[n1;n2;n3;n4]
#Número de grupos
g=size(ng,1)
#Vetor de médias para o grupo 1
MEDIAX1=mean(X1,dims=1)
#Vetor de médias para o grupo 2
MEDIAX2=mean(X2,dims=1)
#Vetor de médias para o grupo 3
MEDIAX3=mean(X3,dims=1)
#Vetor de médias para o grupo 4
MEDIAX4=mean(X4,dims=1)
#Vetor de médias dos grupos
MEDIAS=[MEDIAX1;MEDIAX2;MEDIAX3;MEDIAX4]
#Vetor de médias geral
MEDIAGERAL=mean(X,dims=1)
#Matriz de diferenças das médias dos grupos e a média geral
DIFMEDIA=MEDIAS.-MEDIAGERAL
#Matriz de diferenças das médias dos grupos e a média geral corrigida pelo tamanho dos grupos
DIFMEDIACOR=DIFMEDIA.*ng
178
#Matriz de somas de quadrados dos grupos
B=DIFMEDIACOR'*DIFMEDIA
#Matriz de diferenças dos valores de cada grupo e a média de seu respectivo 
grupo
DIFX1=X1.-MEDIAX1
DIFX2=X2.-MEDIAX2
DIFX3=X3.-MEDIAX3
DIFX4=X4.-MEDIAX4
DIF=[DIFX1;DIFX2;DIFX3;DIFX4]
#Matriz de somas de quadrados dos resíduos
W=DIF'*DIF
#Cálculo do valor de Lambda de Wilks
L=det(W)/det(B+W)
#Cálculo da estatística de Pillai
Pillai=tr(B*inv(B+W))
#Cálculo da estatística de Hotteling
Hotteling=tr(inv(W)*B)
#Cálculo da estatística de Roy
#Autovalores da matriz inversa de W * B
INVWB=inv(W)*B
AV = eigen(INVWB)
#Autovalores ordenados
Valor = reverse(AV.values, dims=1)
ROY=Valor[1,1]/(1+Valor[1,1])
#Comparação do valor de Lambda de Wilks contra o valor crítico de Qui-quadrado
#Cálculo do valor de Lambda de Wilks corrigido
Lcorrigido=-(n-1-((p+g)/2))*log(L)
#Valor crítico de Qui-Quadrado 
using Distributions
QUIQUADRAINV = quantile(Chisq(p * (g - 1)), 0.99)

#MÉTODO 2: MANOVA POR MEIO DE MODELOS DE REGRESSÃO LINEARES MÚLTIPLOS. ALTERNATIVA PARA A INVERSA DE X’X
using DataFrames, CSV, Statistics, LinearAlgebra
#Matriz das variáveis respostas
Y=[
0.50 69.33 27.43 3.24 5113.53
0.53 69.01 27.81 3.18 5115.39
0.51 68.94 27.88 3.19 5043.65
179
0.50 69.07 27.61 3.32 5081.57
0.51 75.24 22.58 2.18 4708.60
0.51 73.24 24.28 2.49 4553.28
0.51 73.89 23.89 2.21 4479.64
0.48 73.93 23.75 2.33 4664.67
0.52 73.44 23.91 2.65 4754.10
0.54 72.15 24.97 2.88 4909.00
0.50 73.81 23.48 2.71 4704.58
0.49 73.80 23.71 2.49 4690.54
0.55 74.63 23.50 1.88 4620.83
0.61 73.63 24.13 2.24 4647.88
0.60 72.24 25.23 2.53 4657.54
0.59 72.69 24.87 2.45 4693.50]
#Matriz de variáveis categóricas para representar os tratamentos
X=[ 
1 1 0 0 0
1 1 0 0 0
1 1 0 0 0
1 1 0 0 0
1 0 1 0 0
1 0 1 0 0
1 0 1 0 0
1 0 1 0 0
1 0 0 1 0
1 0 0 1 0
1 0 0 1 0
1 0 0 1 0
1 0 0 0 1
1 0 0 0 1
1 0 0 0 1
1 0 0 0 1]
#Somas de quadrados e produtos não corrigidas
S=X'X
#Determinação do número de repetições
T=sum(X[:,2:end], dims=1)
T=[0 T]
(m,n)=size(X)
#Matriz nula
Nula=zeros((n-1),n)
#Fator de correção 
A=[
T
Nula]
#Correção da matriz S para admitir inversa
M=S-A
180
#Matriz dos coeficientes
B=inv(M)X'Y
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
#Testes estatísticos
B=SQTrat_Corr
W=SQRes
#Cálculo de Lambda de Wilks
L=det(W)/det(B+W)
#Cálculo da estatística de Pillai
Pillai=tr(B*inv(B+W))
#Cálculo da estatística de Hotteling
Hotteling=tr(inv(W)*B)
#Cálculo da estatística de Roy
#Autovalores da matriz inversa de W * B
INVWB=inv(W)*B
AV = eigen(INVWB)
#Autovalores ordenados
Valor = reverse(AV.values, dims=1)
ROY=Valor[1,1]/(1+Valor[1,1])
