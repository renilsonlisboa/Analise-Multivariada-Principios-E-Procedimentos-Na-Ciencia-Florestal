#Método 2. Análise discriminate para várias populações
#Solução do exemplo 7.2
using DataFrames, CSV, Statistics, LinearAlgebra, Distributions, Plots,
StatsBase;
#Valores das variáveis para cada grupo
X1=[0.47 68.62 28.10 5433.00
0.47 69.02 27.85 5235.00
0.48 69.10 27.67 5092.24
0.51 69.10 27.67 5092.24
0.53 69.13 28.23 5312.00
0.55 69.10 27.67 5092.24
0.45 69.10 27.67 5092.24
225
0.50 69.10 27.67 5092.24
0.54 71.19 25.52 4991.00
0.57 70.68 27.19 5018.00]
X2=[0.49 73.44 24.92 4694.00
0.51 74.27 23.93 4759.00
0.49 74.46 23.53 4870.00
0.57 72.28 24.97 4691.00
0.56 72.37 25.12 4776.00
0.55 70.23 26.77 4878.00
0.52 73.49 23.83 4724.00
0.56 72.39 24.89 4815.00
0.54 73.16 24.37 4705.00
0.51 72.62 23.95 4855.00]
X3=[0.56 80.99 18.50 4434.00
0.58 74.21 23.69 4442.00
0.51 73.77 24.34 4705.00
0.64 76.84 21.61 4613.00
0.54 74.01 23.99 4682.00
0.52 72.19 25.79 4708.00
0.34 71.55 25.49 4683.00
0.59 73.32 24.41 4653.83
0.52 77.53 21.59 4607.00
0.58 74.73 23.30 4607.00]
#Matriz de dados geral
X=[
 X1
 X2
 X3
]
#Número de grupos
g=3
#Matriz de dados padronizados
Z1 = (X1 .- mean(X, dims=1)) ./ std(X, dims=1)
Z2 = (X2 .- mean(X, dims=1)) ./ std(X, dims=1)
Z3 = (X3 .- mean(X, dims=1)) ./ std(X, dims=1)
#Vetores médios dos grupos
MIX1=mean(Z1, dims=1)
MIX2=mean(Z2, dims=1)
MIX3=mean(Z3,dims=1)
#Tamanho das amostras de cada grupo
(n1,p)=size(X1)
(n2,p)=size(X2)
(n3,p)=size(X3)
226
#Vetor médio geral
MI=(n1*MIX1+n2*MIX2+n3*MIX3)/(n1+n2+n3)
#Matrizes de variância-covariância
S1=cov(Z1)
S2=cov(Z2)
S3=cov(Z3)
#Matriz variância-covariância combinada
Sp=((n1-1)*S1+(n2-1)*S2+(n3-1)*S3)/(n1+n2+n3-g)
#Matriz B, representa a variação entre grupos
B=(n1*((MIX1-MI)'*(MIX1-MI)))+(n2*((MIX2-MI)'*(MIX2-MI)))+(n3*((MIX3-
MI)'*(MIX3-MI)))
#Matriz W, representa a variação dentro dos grupos
W=(n1+n2+n3-g)*Sp
#Autovalores e autovetores de WinvB
#Matriz inversa de W e WinvB
WINV=inv(W)
WINVB=WINV*B
#Número de funções discriminantes
s=min(g-1,p)
#Autovalores e autovetores de WINVB
AV = eigen(WINVB)
#Autovalores ordenados
l = reverse(AV.values, dims=1)
l=l[1:s]
#Autovalores acumulados
AVA = cumsum(l)
#Variância explicada acumulada pelas funções discriminantes
VE = AVA ./ AVA[end]*100
#Correlação canônica para cada função discriminante
r=(l./([1].+l)).^0.5
#Teste Lambda de Wilks para cada função discriminante
#Wilks para a primeira função discriminante
Wilks1=(1/(1+l[1,1]))*(1/(1+l[2,1]))
#Wilks para as demais funções discriminantes
Wilks=[1]./([1].+l[2,end])
#Número de grupos
q=g 
227
#Número de observações total
n=n1+n2+n3
#Valores de Qui-quadrado para cada função discriminante
Qui2W1=-(n-1-0.5*(p+q))*log(Wilks1)
Qui2W2=-(n-1-0.5*(p+q))*log.(Wilks)
#Graus de liberdade de cada função discriminante
#m corresponde a ordem do teste de Qui-quadrado
m=1
gl1=(p-m+1)*(q-m)
m=2
gl2=(p-m+1)*(q-m)
#Valores de alfa correspondentes ao Qui-quadrado de cada função discriminante
Alpha1 = ccdf(Chisq(gl1), Qui2W1)
Alpha2 = ccdf(Chisq(gl2), Qui2W2)
#Obtenção dos coeficientes das funções discriminantes 
#Matriz B0, sem n
B0=(((MIX1-MI)'*(MIX1-MI)))+(((MIX2-MI)'*(MIX2-MI)))+(((MIX3-MI)'*(MIX3-MI)))
#Matriz inversa de W multiplicada por B0
WINVB0=WINV*B0
#Autovalores e autovetores de WINVB0
AV = eigen(WINVB0)
#Autovalores ordenados
l = reverse(AV.values, dims=1)
#Autovetores correspondentes aos autovalores
e = reverse(AV.vectors, dims=2)
#Fator de escalonamento e autovetores escalonados
e1=e[:,1]
e2=e[:,2]
k1=1/sqrt((e1)'*Sp*e1)
e1escalonado=k1*e1
k2=1/sqrt((e2)'*Sp*e2)
e2escalonado=k2*e2
#Coeficientes das funções discriminantes não padronizadas
ee=[e1escalonado e2escalonado]
#Coeficientes das funções discriminantes padronizados
Sp
ee
D=diag(Sp)
DR=sqrt.(D)
228
eep=DR.*ee
#Matriz de estruturas
Spd=diagm(DR)
Rp=inv(Spd)*Sp*inv(Spd) #Matriz de correlação
#Matriz de estruturas para a primeira função discriminante
ee1=ee[:,1]
c1=ee1'*Spd
me1=Rp*c1'
#Matriz de estruturas para a segunda função discriminante
ee2=ee[:,2]
c2=ee2'*Spd
me2=Rp*c2'
#Matriz de estruturas final
me=[me1 me2] 
#Gráfico das funções discriminantes padronizadas
#Grupo 1
G1y1= Z1 * eep[:,1]
G1y2= Z1 * eep[:,2]
scatter(G1y1, G1y2, label="Grupo 1")
#Grupo 2
G2y1= Z2 * eep[:,1]
G2y2= Z2 * eep[:,2]
scatter!(G2y1, G2y2, label="Grupo 2")
#Grupo 3
G3y1= Z3 * eep[:,1]
G3y2= Z3 * eep[:,2]
scatter!(G3y1, G3y2, label="Grupo 3")
#Centróide grupo 1
CG1y1=mean(G1y1)
CG1y2=mean(G1y2)
CG1=[CG1y1 CG1y2]
#Centróide grupo 2
CG2y1=mean(G2y1)
CG2y2=mean(G2y2)
CG2=[CG2y1 CG2y2]
#Centróide grupo 3
CG3y1=mean(G3y1)
CG3y2=mean(G3y2)
CG3=[CG3y1 CG3y2]
229
#Solução do exemplo 7.3
#Obtencão das informações para a matriz de confusão
#Distância euclidiana do grupo 1 em relação ao centróide do grupo 1
dG1C1=sqrt.((G1y1.-CG1y1).^2+(G1y2.-CG1y2).^2)
#Distância euclidiana do grupo 1 em relação ao centróide do grupo 2
dG1C2=sqrt.((G1y1.-CG2y1).^2+(G1y2.-CG2y2).^2)
#Distância euclidiana do grupo 1 em relação ao centróide do grupo 3
dG1C3=sqrt.((G1y1.-CG3y1).^2+(G1y2.-CG3y2).^2)
#matriz de distância do grupo 1 em relação aos centróides
DG1=[dG1C1 dG1C2 dG1C3]
# Índice do centróide mais próximo
col_min_G1 = mapslices(DG1 -> argmin(DG1), DG1; dims = 2)
# Converter para vetor
GrupoClassificado1 = vec(col_min_G1)
# Contagem de ocorrências para o grupo 1
contagem1 = countmap(col_min_G1)
#Distância euclidiana do grupo 2 em relação ao centróide do grupo 1
dG2C1=sqrt.((G2y1.-CG1y1).^2+(G2y2.-CG1y2).^2)
#Distância euclidiana do grupo 2 em relação ao centróide do grupo 2
dG2C2=sqrt.((G2y1.-CG2y1).^2+(G2y2.-CG2y2).^2)
#Distância euclidiana do grupo 2 em relação ao centróide do grupo 3
dG2C3=sqrt.((G2y1.-CG3y1).^2+(G2y2.-CG3y2).^2)
#matriz de distância do grupo 2 em relação aos centróides
DG2=[dG2C1 dG2C2 dG2C3]
#Índice do centróide mais próximo
col_min_G2 = mapslices(DG2 -> argmin(DG2), DG2; dims = 2)
# Converter para vetor
GrupoClassificado2 = vec(col_min_G2)
# Contagem de ocorrências para o grupo 2
contagem2 = countmap(col_min_G2)
#Distância euclidiana do grupo 3 em relação ao centróide do grupo 1
dG3C1=sqrt.((G3y1.-CG1y1).^2+(G3y2.-CG1y2).^2)
#Distância euclidiana do grupo 3 em relação ao centróide do grupo 2
dG3C2=sqrt.((G3y1.-CG2y1).^2+(G3y2.-CG2y2).^2)
#Distância euclidiana do grupo 3 em relação ao centróide do grupo 3
dG3C3=sqrt.((G3y1.-CG3y1).^2+(G3y2.-CG3y2).^2)
#matriz de distância do grupo 3 em relação aos centróides
DG3=[dG3C1 dG3C2 dG3C3]
#Índice do centróide mais próximo
col_min_G3 = mapslices(DG3 -> argmin(DG3), DG3; dims = 2)
# Converter para vetor
GrupoClassificado3 = vec(col_min_G3)
# Contagem de ocorrências para o grupo 3
contagem3 = countmap(col_min_G3)
230
#Classificação de uma nova observação
Xnovo=[0.62 76.78 21.62 4614.11]
#Dados padronizados da nova observação
Znovo = (Xnovo .- mean(X, dims=1)) ./ std(X, dims=1)
#Escores das funções discriminantes padronizadas
Ny1= Znovo * eep[:,1]
Ny2= Znovo * eep[:,2]
#Distância euclidiana da nova observação em relação ao centróide do grupo 1
dNovoC1=sqrt.((Ny1.-CG1y1).^2+(Ny2.-CG1y2).^2)
#Distância euclidiana da nova observação em relação ao centróide do grupo 2
dNovoC2=sqrt.((Ny1.-CG2y1).^2+(Ny2.-CG2y2).^2)
#Distância euclidiana da nova observação em relação ao centróide do grupo 3
dNovoC3=sqrt.((Ny1.-CG3y1).^2+(Ny2.-CG3y2).^2)
#matriz de distância da nova observação em relação aos centróides
DNovo=[dNovoC1 dNovoC2 dNovoC3]
# Índice do centróide mais próximo
col_min_novo = mapslices(DNovo -> argmin(DNovo), DNovo; dims = 2)
# Converter para vetor
GrupoClassificadonovo = vec(col_min_novo)