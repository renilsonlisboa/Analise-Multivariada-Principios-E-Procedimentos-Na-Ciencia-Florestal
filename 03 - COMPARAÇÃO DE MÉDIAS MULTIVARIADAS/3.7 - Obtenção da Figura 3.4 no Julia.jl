# OBTENÇÃO DA FIGURA 3.4 NO JULIA
# GRÁFICO DA DISTRIBUIÇÃO NORMAL BIVARIADA PARA AS VARIÁVEIS 
# DIÂMETRO À ALTURA DO PEITO E ALTURA COM CORRELAÇÃO ALTA (r = 0.88)
using LinearAlgebra, Distributions, PyPlot
# VETOR DE MÉDIAS CENTRALIZADO
mu = [0.0, 0.0]
# MATRIZ VARIÂNCIA-COVARIÂNCIA DAS VARIÁVEIS ORIGINAIS
S = [13.38 9.80;
 9.80 9.31]
# ORDENADAS PARA AS VARIÁVEIS
x1 = -12:0.2:12
x2 = -12:0.2:12
# PLOTAGEM DO GRID (equivalente ao meshgrid do MATLAB)
X1 = repeat(collect(x1)', length(x2), 1)
X2 = repeat(collect(x2), 1, length(x1))
# ORGANIZAÇÃO DOS PONTOS EM MATRIZ (equivalente a X=[X1(:) X2(:)])
X = hcat(vec(X1), vec(X2))
# CÁLCULO DOS VALORES DE DENSIDADE PROBABILÍSTICA
# (equivalente ao mvnpdf do MATLAB)
dist = MvNormal(mu, S)
y = pdf.(dist, eachrow(X))
# REORGANIZAÇÃO PARA MATRIZ (reshape como no MATLAB)
y = reshape(y, length(x2), length(x1))
# OBTENÇÃO DO GRÁFICO
figure()
surf(x1, x2, y)
axis([-12, 12, -12, 12, 0, 0.03])
xlabel("DAP")
ylabel("ALTURA")
zlabel("FDP")
# GRÁFICO DA DISTRIBUIÇÃO NORMAL BIVARIADA PARA AS VARIÁVEIS 
# TEOR DE CARBONO NO FUSTE E ALTURA DO FUSTE COM CORRELAÇÃO BAIXA (r = -0,056)
using LinearAlgebra, Distributions, PyPlot
# VETOR DE MÉDIAS CENTRALIZADO
mu = [0.0, 0.0]
# MATRIZ VARIÂNCIA-COVARIÂNCIA DAS VARIÁVEIS ORIGINAIS
S = [0.035 -0.032;
 -0.032 3.312]
# ORDENADAS PARA AS VARIÁVEIS
x1 = -0.7:0.02:0.7
x2 = -6:0.2:6
# PLOTAGEM DO GRID (equivalente ao meshgrid do MATLAB)
X1 = repeat(collect(x1)', length(x2), 1)
X2 = repeat(collect(x2), 1, length(x1))
# ORGANIZAÇÃO DOS PONTOS EM MATRIZ (equivalente a X=[X1(:) X2(:)])
X = hcat(vec(X1), vec(X2))
# CÁLCULO DOS VALORES DE DENSIDADE PROBABILÍSTICA
# (equivalente ao mvnpdf do MATLAB)
dist = MvNormal(mu, S)
y = pdf.(dist, eachrow(X))
68
# REORGANIZAÇÃO PARA MATRIZ (reshape como no MATLAB)
y = reshape(y, length(x2), length(x1))
# OBTENÇÃO DO GRÁFICO
figure()
surf(x1, x2, y)
axis([-0.7, 0.7, -6, 6, 0, 0.5])
xlabel("TEOR DE C NO FUSTE")
ylabel("ALTURA DO FUSTE")
zlabel("PDF")