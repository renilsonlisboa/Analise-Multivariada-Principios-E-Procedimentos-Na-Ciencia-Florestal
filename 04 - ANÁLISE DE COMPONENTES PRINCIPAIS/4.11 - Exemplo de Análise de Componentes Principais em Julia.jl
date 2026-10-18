using DataFrames, CSV
#Dados originais
X=CSV.read("Z:\\Alexandre\\Julia\\Multivariada\\dendrometria.csv", 
DataFrame)
#Matriz X com dados originais 
X=[X.d X.h X.g X.f];

#Componentes principais usando matriz variância-covariância
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
#Variância explicada acumulada pelas componentes principais
VE = AVA ./ AVA[end]
#Cálculo dos loadings
#a)Raiz quadrada dos autovalores
RaizValor=Valor.^0.5
#b)Desvio padrão da matriz variância-covariância
Desvio=diag(S).^0.5
#c)Número de variáveis
n=length(Desvio)
#d)Correlação das componentes principais com as variáveis originais -
loadings
Loadings = [Vetor[c,l]*RaizValor[l]/Desvio[c] for l in 1:n, c in 1:n]'
#Escores das componentes principais 
Escores = X * Vetor
#Componentes principais
CP1=Escores[:,1];
122
CP2=Escores[:,2];
CP3=Escores[:,3];
CP4=Escores[:,4];
#Organização dos scores das duas primeiras componentes principais em uma 
tabela
using JuliaDB
t = table((pc1 = CP1, pc2 = CP2))
#Gráfico entre as duas primeiras componentes principais
using StatsPlots
scatter(CP1,CP2)
xlabel!("CP1")
ylabel!("CP2")
#Componentes principais usando a matriz de correlação (dados padronizados - Z)
#Matriz de dados padronizados
Z = (X .- mean(X, dims=1)) ./ std(X, dims=1)
#Matriz de correlação dos dados originais
CORRX=cor(X)
#Matriz de correlação dos dados padronizados
CORRNORMAL=cor(Z)
#Matriz variância-covariância dos dados padronizados
COVS=cov(Z)
#Autovalores e autovetores da matriz de correlação
AVN = eigen(CORRX)
#Autovalores ordenados
ValorN = reverse(AVN.values, dims=1)
#Autovetores ordenados por autovalores
VetorN = reverse(AVN.vectors, dims=2)
#Autovalores acumulados
AVAN = cumsum(ValorN)
#Variância explicada acumulada pelas componentes principais
VEN = AVAN ./ tr(CORRX)
#Cálculo dos loadings
123
#a)Raiz quadrada dos autovalores
RaizValor=ValorN.^0.5
#b)Número de variáveis
n=length(RaizValor)
#c)Correlação das componentes principais com as variáveis originais -
loadings
LoadingsN = [VetorN[c,l]*RaizValor[l] for l in 1:n, c in 1:n]'
# Escores das componentes principais padronizadas
EscoresN = Z * VetorN
#Componentes principais
CP1N=EscoresN[:,1];
CP2N=EscoresN[:,2];
CP3N=EscoresN[:,3];
CP4N=EscoresN[:,4];
#Gráfico entre as duas primeiras componentes principais padronizadas
using StatsPlots
scatter(CP1N,CP2N)
xlabel!("CP1")
ylabel!("CP2")
#Gráfico do teste Scree 
plot(ValorN, legend = false)
xlabel!("Componente Principal")
ylabel!("Autovalor")
# Importância das variáveis na primeira componente principal padronizada
plot(LoadingsN[:,1], seriestype = :bar, legend = false)
ylims!((-1,1))
xlabel!("Variável")
ylabel!("Correlação")
title!("Primeira componente principal")
# Importância das variáveis em cada componente principal padronizada
groupedbar(Loadings, bar_position = :dodge, bar_width=0.7, legend=:top)
ylims!((-1,1))
xlabel!("Variável")
ylabel!("Correlação")