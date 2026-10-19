using LinearAlgebra; Statistics
# Dados de entrada 
X = [
 -1.0200 -1.3000 1.3000 2.3600
 -0.9500 -1.1600 1.1900 1.5700
 -0.8900 -1.1300 1.1100 0.9900
 -0.3100 -1.1300 1.1100 0.9900
 0.1300 -1.1200 1.3600 1.8800
 0.5300 -1.1300 1.1100 0.9900
 -1.3500 -1.1300 1.1100 0.9900
 -0.4700 -1.1300 1.1100 0.9900
 0.2700 -0.4200 0.1700 0.5900
 0.8200 -0.5900 0.9000 0.7000
 -0.7300 0.3600 -0.1000 -0.6100
 -0.2000 0.6400 -0.5300 -0.3500
 -0.6100 0.7100 -0.7100 0.1000
 0.7600 -0.0400 -0.0700 -0.6200
 0.6100 -0.0100 -0.0100 -0.2800
 0.4600 -0.7500 0.7200 0.1300
 -0.1700 0.3800 -0.5800 -0.4900
 0.7500 0 -0.1100 -0.1200
 0.3600 0.2600 -0.3400 -0.5600
 -0.2100 0.0800 -0.5200 0.0400
 0.7000 2.9500 -2.9200 -1.6500
 0.9500 0.6200 -0.6400 -1.6200
 -0.2900 0.4700 -0.3500 -0.5600
 2.1900 1.5300 -1.5500 -0.9300
 0.2400 0.5500 -0.5100 -0.6600
 -0.1400 -0.0700 0.2800 -0.5500
 -3.3900 -0.2900 0.1500 -0.6500
253
 1.1800 0.3200 -0.3200 -0.7700
 -0.1700 1.7600 -1.5600 -0.9600
 0.9700 0.8000 -0.8100 -0.9600
]
n, p = size(X)
# Definição do número de grupos
k = 3
# Inicialização dos centroides
# Índices escolhidos arbitrariamente
centroides = X[[1, 11, 21], :] 
# Função distância euclidiana
function dist_euclidiana(x, c)
 return sqrt(sum((x .- c).^2))
end
# Iterações do k-means
maxiter = 100
tol = 1e-6
grupos = zeros(Int, n)
centroides_antigos = copy(centroides)
for iter in 1:maxiter
 # PASSO 1: Atribuição aos grupos
 for i in 1:n
 distancias = zeros(k)
 for g in 1:k
 distancias[g] = dist_euclidiana(X[i, :], centroides[g, :])
 end
 grupos[i] = argmin(distancias)
 end
 # PASSO 2: Atualização dos centroides
 for g in 1:k
 idx = findall(grupos .== g)
 centroides[g, :] = mean(X[idx, :], dims=1)
 end
 # PASSO 3: Critério de parada
 deslocamento = maximum(abs.(centroides .- centroides_antigos))
 if deslocamento < tol
 println("Convergência atingida na iteração ", iter)
 break
 end
 centroides_antigos = copy(centroides)
end
#Resultados finais
println("\nGrupos atribuídos a cada observação:")
println(grupos)
println("\nCentroides finais:")
254
println(centroides)
# Número de observações em cada grupos
println("Número de observações em cada grupo:")
for g in 1:k
 println("Grupo $g: ", sum(grupos .== g), " observações")
end