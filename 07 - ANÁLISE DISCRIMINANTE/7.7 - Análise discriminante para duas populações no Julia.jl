#Método 1. Análise discriminante para duas populações
#Solução do exemplo 7.1
using DataFrames, CSV, Statistics, LinearAlgebra, Distributions, Plots;
#Matriz de dados da população 1
X1=[0.47 68.62 28.10 3.29 5433.00
0.47 69.02 27.85 3.14 5235.00
0.48 69.10 27.67 3.23 5092.24
0.51 69.10 27.67 3.23 5092.24
0.53 69.13 28.23 2.64 5312.00
223
0.55 69.10 27.67 3.23 5092.24
0.45 69.10 27.67 3.23 5092.24
0.50 69.10 27.67 3.23 5092.24
0.54 71.19 25.52 3.29 4991.00
0.57 70.68 27.19 2.13 5018.00
0.50 69.17 27.44 3.39 5157.00
0.51 69.10 27.67 3.23 5092.24
0.66 69.10 27.67 3.23 5092.24
0.50 67.48 29.47 3.05 5121.00
0.44 69.10 27.67 3.23 5092.24
0.48 68.63 27.51 3.86 5168.00
0.58 69.10 27.67 3.23 5092.24
0.40 69.75 26.33 3.92 5112.00
0.41 70.67 25.99 3.34 4945.00
0.42 70.30 26.14 3.57 5059.00
0.56 69.41 27.19 3.40 5003.00]
#Matriz de dados da população 2
X2=[0.40 75.98 23.08 0.94 4563.00
0.44 77.52 21.45 1.04 4601.00
0.40 75.94 22.69 1.37 4568.00
0.43 75.21 23.40 1.40 4518.57
0.45 75.97 23.14 0.89 4659.00
0.45 75.85 22.89 1.25 4552.00
0.46 75.42 23.44 1.14 4508.00
0.48 77.08 21.93 0.99 4481.00
0.44 77.53 21.37 1.10 4522.00
0.41 77.07 22.05 0.88 4412.00
0.42 76.05 22.47 1.47 4449.00
0.43 75.21 23.40 1.40 4518.57
0.44 75.07 23.59 1.35 4549.00
0.57 75.21 23.40 1.40 4620.00
0.32 73.64 24.73 1.63 4541.00
0.35 75.02 23.58 1.40 4469.00
0.44 75.66 22.92 1.42 4411.00
0.39 75.64 23.38 0.98 4487.00
0.34 75.20 23.32 1.49 4406.00
0.47 77.08 21.88 1.03 4386.00
0.40 78.15 20.71 1.14 4270.00]
#Matriz de dados geral
X=[X1
X2]
# Médias das populações
MI1 = mean(X1,dims=1)
MI2 = mean(X2,dims=1)
# Matriz variância-covariância
S = cov(X)
224
# Função Discriminante Linear de Fisher
Y = (MI1 - MI2) * inv(S)
# Médias das populações univariadas
MI1Y = Y * MI1'
MI2Y = Y * MI2'
# Variância das populações univariadas
VY = Y * S * Y'
# Ponto médio entre as populações univariadas
m = (MI1Y + MI2Y) / 2
#Gráfico das funções univariadas
x = -2130:0.01:-2090
a=MI1Y[1,1]
b=MI2Y[1,1]
v=VY[1,1]
y1 = (1 / (sqrt(v) * sqrt(2*pi))) .*
 exp.(-0.5 .* (((x .- (a)) ./ sqrt(v)).^2))
y2 = (1 / (sqrt(v) * sqrt(2*pi))) .*
 exp.(-0.5 .* (((x .- b) ./ sqrt(v)).^2))
plot(x, y1, label="População 1")
plot!(x, y2, label="População 2")
#Probabilidade da classificação errada
z=(m- MI1Y)/sqrt(v)
AREA=cdf(Normal(0,1),z)
p=AREA*2
#Alocação de uma nova observação
X0=[0.45 68.3 26.4 3.2 5100.41]
Y0=Y*X0'