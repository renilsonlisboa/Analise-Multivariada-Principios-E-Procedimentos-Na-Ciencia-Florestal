#Análise de imagem
using Statistics, LinearAlgebra
#Padronização dos dados originais
Z = (X .- mean(X, dims=1)) ./ std(X, dims=1)
#Matriz de correlação dos dados originais
R=cor(X)
#Inversa da matriz de correlação R
RINV=inv(R)
#Matriz recíproca da inversa RINV
REC=RINV.^-1
#Diagonal da recíproca da inversa
D=diag(REC)
D=diagm(D)
B=RINV*D
#Matriz identidade
II=Matrix(I, 4, 4)
#Matriz modificada de pesos - matriz A
A=II-B
#Matriz comum
CC=Z*A
#Matriz variância-covariância da matriz comum
VARCOV_CC=cov(CC)
#Matriz correlação da matriz comum
R_CC=cor(CC)
#Autovalores e autovetores da matriz de correlação da matriz comum
ValVet=eigen(R_CC)
#Autovalores ordenados
Valor = reverse(ValVet.values, dims=1)
160
#Autovetores ordenados por autovalores - matriz modal
Modal = reverse(ValVet.vectors, dims=2)
#Autovalores acumulados
AVA = cumsum(Valor)
#Variância explicada acumulada pelos fatores
VE = AVA ./ AVA[end]
#Matriz de dados transformados 
Y = CC * Modal
#Matriz variância-covariância da matriz de dados transformados - Y
VarCov_Y=cov(Y)
#Matrix de correlação da matriz de dados transformados - Y
R_Y=cor(Y)
#Raíz quadrada dos autovalores
Raiz_Valor=Valor.^0.5
#Cáluclo dos carregamentos - loadings
n=length(R[:,1])
L = [Modal[c,l]*Raiz_Valor[l] for l in 1:n, c in 1:n]'
