# ANÁLISE DE AGRUPAMENTO COM DISTÂNCIA EUCLIDIANA E LIGAÇÃO SIMPLES
using LinearAlgebra; Statistics
# DADOS DE ENTRADA
X = [
 -1.0200 -1.3000 1.3000 2.3600;
 -0.9500 -1.1600 1.1900 1.5700;
 -0.8900 -1.1300 1.1100 0.9900;
 -0.3100 -1.1300 1.1100 0.9900;
 0.1300 -1.1200 1.3600 1.8800;
 0.5300 -1.1300 1.1100 0.9900;
 -1.3500 -1.1300 1.1100 0.9900;
 -0.4700 -1.1300 1.1100 0.9900;
 0.2700 -0.4200 0.1700 0.5900;
 0.8200 -0.5900 0.9000 0.7000;
 -0.7300 0.3600 -0.1000 -0.6100;
 -0.2000 0.6400 -0.5300 -0.3500;
 -0.6100 0.7100 -0.7100 0.1000;
 0.7600 -0.0400 -0.0700 -0.6200;
 0.6100 -0.0100 -0.0100 -0.2800;
 0.4600 -0.7500 0.7200 0.1300;
 -0.1700 0.3800 -0.5800 -0.4900;
 0.7500 0.0000 -0.1100 -0.1200;
 0.3600 0.2600 -0.3400 -0.5600;
 -0.2100 0.0800 -0.5200 0.0400;
 0.7000 2.9500 -2.9200 -1.6500;
250
 0.9500 0.6200 -0.6400 -1.6200;
 -0.2900 0.4700 -0.3500 -0.5600;
 2.1900 1.5300 -1.5500 -0.9300;
 0.2400 0.5500 -0.5100 -0.6600;
 -0.1400 -0.0700 0.2800 -0.5500;
 -3.3900 -0.2900 0.1500 -0.6500;
 1.1800 0.3200 -0.3200 -0.7700;
 -0.1700 1.7600 -1.5600 -0.9600;
 0.9700 0.8000 -0.8100 -0.9600
]
numero_observacoes = size(X, 1)
# MATRIZ DE DISTÂNCIAS EUCLIDIANAS
matriz_distancias = zeros(numero_observacoes, numero_observacoes)
for i in 1:numero_observacoes
 for j in i+1:numero_observacoes
 matriz_distancias[i,j] = norm(X[i,:] - X[j,:])
 matriz_distancias[j,i] = matriz_distancias[i,j]
 end
end
# INICIALIZAÇÃO
agrupamentos = [[i] for i in 1:numero_observacoes]
distancias_atuais = copy(matriz_distancias)
matriz_ligacao = zeros(numero_observacoes - 1, 3)
matriz_cofenetica = zeros(numero_observacoes, numero_observacoes)
etapa = 0
# ALGORITMO HIERÁRQUICO (LIGAÇÃO SIMPLES)
while length(agrupamentos) > 1
 etapa += 1
 # IGNORAR A DIAGONAL
 matriz_temporaria = copy(distancias_atuais)
 for k in 1:size(matriz_temporaria, 1)
 matriz_temporaria[k,k] = Inf
 end
 # MENOR DISTÂNCIA
 indice = argmin(matriz_temporaria)
 i, j = Tuple(indice)
 menor_distancia = matriz_temporaria[i,j]
 agrupamento_A = agrupamentos[i]
 agrupamento_B = agrupamentos[j]
 # ATUALIZAR MATRIZ COFENÉTICA
 for p in agrupamento_A
 for q in agrupamento_B
 if matriz_cofenetica[p,q] == 0.0
 matriz_cofenetica[p,q] = menor_distancia
 matriz_cofenetica[q,p] = menor_distancia
 end
 end
 end
251
 # MATRIZ DE LIGAÇÃO
 matriz_ligacao[etapa, :] = [
 agrupamento_A[1],
 agrupamento_B[1],
 menor_distancia]
 # NOVO AGRUPAMENTO
 novo_agrupamento = vcat(agrupamento_A, agrupamento_B)
 # REMOVER AGRUPAMENTOS ANTIGOS
 # (ORDEM DECRESCENTE)
 for indice_agrupamento in sort([i, j], rev=true)
 deleteat!(agrupamentos, indice_agrupamento)
 distancias_atuais = distancias_atuais[
 setdiff(1:end, indice_agrupamento),
 setdiff(1:end, indice_agrupamento)]
 end
 # CALCULAR NOVAS DISTÂNCIAS
 # (LIGAÇÃO SIMPLES)
 numero_agrupamentos = length(agrupamentos)
 novas_distancias = zeros(numero_agrupamentos)
 for k in 1:numero_agrupamentos
 lista_distancias = Float64[]
 for a in novo_agrupamento
 for b in agrupamentos[k]
 push!(
 lista_distancias,
 matriz_distancias[a,b])
 end
 end
 novas_distancias[k] = minimum(lista_distancias)
 end
 # ATUALIZAR MATRIZ DE DISTÂNCIAS
 distancias_atuais = [
 distancias_atuais novas_distancias;
 novas_distancias' 0.0]
 push!(agrupamentos, novo_agrupamento)
end
# CORRELAÇÃO COFENÉTICA
mascara = triu(trues(numero_observacoes, numero_observacoes), 1)
distancias_originais = matriz_distancias[mascara]
distancias_cofeneticas = matriz_cofenetica[mascara]
correlacao_cofenetica = cor(distancias_originais, distancias_cofeneticas)
println("Correlação cofenética = ", correlacao_cofenetica)
252
# MAIOR DISTÂNCIA NO DENDROGRAMA
maior_distancia = maximum(matriz_ligacao[:,3])
println("Maior distância no dendrograma = ", maior_distancia)
# OBTENÇÃO DO DENDROGRAMA
using Clustering; Distances; StatsBase; StatsPlots; LinearAlgebra; Statistics
# MATRIZ DE DISTÂNCIAS EUCLIDIANAS
matriz_distancias = pairwise(Euclidean(), X; dims=)
# AGRUPAMENTO HIERÁRQUICO
# MÉTODO DE LIGAÇÃO SIMPLES
agrupamento_hierarquico = hclust(matriz_distancias, linkage=:single)
# REPRESENTAÇÃO GRÁFICA DO DENDROGRAMA
plot(agrupamento_hierarquico)