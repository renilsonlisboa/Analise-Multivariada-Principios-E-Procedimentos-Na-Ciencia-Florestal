using DataFrames, CSV
#Dados originais
X=CSV.read("Z:\\Alexandre\\Julia\\Multivariada\\dendrometria.csv", 
DataFrame)
155
#Matriz X com dados originais 
X=[X.d X.h X.g X.f];
#Análise de fatores usando matriz variância-covariância - dados originais
using Statistics, LinearAlgebra
#Vetor de médias
Média=mean(X,dims=1)
#Vetor de desvios-padrão
Desvio=std(X,dims=1)
#Matriz variância-covariância
S=cov(X)
#Variância total
TraçoS=tr(S)
#Autovalores e autovetores da matriz variância-covariância
AV = eigen(S)
#Autovalores ordenados
Valor = reverse(AV.values, dims=1)
#Autovetores ordenados por autovalores
Vetor = reverse(AV.vectors, dims=2)
#Autovalores acumulados
AVA = cumsum(Valor)
#Variância explicada acumulada pelos fatores
VE = AVA ./ AVA[end]
#Cálculo dos loadings
#a)Raiz quadrada dos autovalores
RaizValor=Valor.^0.5
#b)Desvio padrão da matriz variância-covariância
Desvio=diag(S).^0.5
#c)Número de variáveis
n=length(Desvio)
#d)Correlação dos fatores com as variáveis originais - loadings
Loadings = [Vetor[c,l]*RaizValor[l]/Desvio[c] for l in 1:n, c in 1:n]'
#Escolha do número de fatores - critério -> autovalores maiores do que um
Valores=filter(x -> x > 1, Valor) 
156
#Número de fatores
NF=length(Valores)
#Vetores associados aos fatores selecionados
Vetores=Vetor[:,1:NF]
#Loadings associados aos fatores selecionados
L=Loadings[:,1:NF]
#Comunalidades
L2=L.^2
Comunalidades= sum(L2, dims=2)
#Variâncias específicas
s=diag(S)
VEspecifica=s-Comunalidades
#Matriz de resíduos
MatrizResíduos=S-(L*L'+diagm(VEspecifica[:]))
#Escores para os fatores selecionados
FE=X*L
#Análise de fatores usando a matriz de correlação - dados padronizados Z
#Matriz de dados padronizados
Z = (X .- mean(X, dims=1)) ./ std(X, dims=1)
#Matriz de correlação dos dados originais
CORRX=cor(X)
#Autovalores e autovetores da matriz de correlação
AVN = eigen(CORRX)
#Autovalores ordenados
ValorN = reverse(AVN.values, dims=1)
#Autovetores ordenados por autovalores
VetorN = reverse(AVN.vectors, dims=2)
#Autovalores acumulados
AVAN = cumsum(ValorN)
#Variância explicada acumulada pelos fatores
VEN = AVAN ./ tr(CORRX)
#Cálculo dos loadings
#a)Raiz quadrada dos autovalores
RaizValorN=ValorN.^0.5
157
#b)Número de variáveis
n=length(RaizValorN)
#c)Correlação dos fatores com as variáveis originais - loadings
LoadingsN = [VetorN[c,l]*RaizValorN[l] for l in 1:n, c in 1:n]'
#Escolha do número de fatores - critério -> número de fatores que explicam até 80% da variância total
ValoresN=filter(x -> x < 0.8, VEN)
#Número de fatores
NFN=length(ValoresN)
#Vetores associados aos fatores selecionados
VetoresN=VetorN[:,1:NFN]
#Loadings associados aos fatores selecionados
LN=LoadingsN[:,1:NFN]
#Comunalidades
L2N=LN.^2
ComunalidadesN= sum(L2N, dims=2)
#Variâncias específicas
r=diag(CORRX)
VEspecificaN=r-ComunalidadesN
#Matriz de resíduos
MatrizResíduosN=CORRX-(LN*LN'+diagm(VEspecificaN[:]))
#Escores padronizados para os fatores selecionados
FN = Z * LN
#Gráfico dos escores padronizados entre os fatores selecionados
using StatsPlots
scatter(FN[:,1], FN[:,2], legend = false)
ylims!((-4,4))
xlims!((-4,4))
xlabel!("Fator 1")
ylabel!("Fator 2")
#Gráfico dos fatores selecionados
scatter(LN[:,1], LN[:,2], framestyle =:zerolines,legend = false)
ylims!((-1,1))
xlims!((-1,1))
xlabel!("Fator 1")
ylabel!("Fator 2")
158
#Rotação dos fatores
#ROTAÇÃO VARIMAX
using LinearAlgebra
X=[0.7814 0.0862
0.8932 -0.06777
0.4335 0.8215
-0.6642 0.5465]
#Passo 1 – Matriz de loadings
X # matriz de loadings não rotacionados
#Passo 2 – Dimensões da matriz de loadings
NV, NF = size(X)
#Passo 3 – Matriz transposta de loadings
XT = transpose(X)
#Passo 4 – Inicialização da matriz de rotação (T)
T = Matrix{Float64}(I, NF, NF)
#Passo 5 – Definição do número de iterações
#Número máximo de iterações
NMI = 200
ITER = 0
#Passo 6 – Início do processo iterativo
while ITER < NMI
#Passo 7 – Loadings atuais rotacionados
LR = X * T
#Passo 8 – Potência cúbica (maximiza a variância dos quadrados dos loadings)
PC = LR .^ 3
#Passo 9 – Matriz diagonal F
F = (1 / NV) * Diagonal(diag(transpose(LR) * LR))
#Passo 10 – Matriz GR (gradiente da rotação)
#Direção de ajuste da matriz de rotação para maximizar o critério Varimax.
GR = XT * (PC - LR * F)
#Passo 11 – Atualização da matriz de rotação
T = inv(sqrt(GR * transpose(GR))) * GR
#Passo 12 – Atualização do contador
ITER += 1
end
#Passo 13 – Loadings finais rotacionados
LR = X * T 
159
#Gráfico para representar os fatores e os fatores rotacionados
scatter(LN[:,1], LN[:,2], framestyle =:zerolines, label = "L")
scatter!(LR[:,1], LR[:,2], framestyle =:zerolines, label = "L*")
ylims!((-1,1))
xlims!((-1,1))
xlabel!("Fator 1")
ylabel!("Fator 2")