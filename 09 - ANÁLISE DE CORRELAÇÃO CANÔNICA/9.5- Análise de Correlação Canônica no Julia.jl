# Análise de correlação canônica
# Pacotes necessários
using LinearAlgebra; Statistics; Distributions 
# Entrada dos dados
# Matriz de dados X (n × (p+q))
X = [
 174.45 31.64 206.09 19.74 20.20 7.70 0.65
84.41 8.76 93.17 15.12 19.10 4.90 0.50
60.97 4.91 65.88 12.89 17.30 3.00 0.57
47.25 2.06 49.30 11.62 16.90 7.20 0.52
114.06 27.95 142.00 16.55 20.50 11.80 0.56
157.39 37.48 194.87 19.10 21.90 7.00 0.60
72.03 15.14 87.17 14.32 18.40 9.20 0.56
88.98 15.66 104.63 15.12 19.30 12.10 0.58
45.29 5.96 51.25 11.94 16.70 9.50 0.50
101.20 18.46 119.65 15.28 19.50 7.70 0.59
94.16 12.22 106.38 16.23 18.30 4.60 0.52
83.79 10.13 93.92 15.12 19.30 3.80 0.59
200.00 54.52 254.51 20.53 19.50 9.50 0.58
49.66 6.04 55.70 13.21 18.40 2.80 0.57
20.96 3.59 24.55 8.75 14.40 11.20 0.48
71.90 15.34 87.24 13.05 18.00 3.50 0.62
115.57 21.25 136.82 16.23 19.60 9.00 0.62
69.30 11.50 80.80 12.25 17.90 3.80 0.63
16.46 2.49 18.94 7.64 9.90 6.70 0.51
107.51 19.67 127.17 15.44 19.20 6.80 0.56
]
# Definição dos conjuntos de variáveis
# X1: primeiro conjunto (n × p)
X1 = X[:, 1:3]
# X2: segundo conjunto (n × q)
X2 = X[:, 4:7]
# Número de observações
n = size(X1, 1)
# Matrizes de variância–covariância
VC11 = cov(X1) # Covariância de X1 (p × p)
VC22 = cov(X2) # Covariância de X2 (q × q)
VC12 = cov(X1, X2) # Covariância cruzada (p × q)
VC21 = VC12' # Transposta (q × p)
VC = cov(X) # Covariância total
278
# Determinação dos coeficientes canônicos
# Matriz associada aos coeficientes a
PARESA = inv(sqrt(VC11)) * VC12 * inv(VC22) * VC21 * inv(sqrt(VC11))
# Autovalores e autovetores de a
EigA = eigen(PARESA)
# Coeficientes canônicos a
a = inv(sqrt(VC11)) * EigA.vectors
# Matriz associada aos coeficientes b
PARESB = inv(sqrt(VC22)) * VC21 * inv(VC11) * VC12 * inv(sqrt(VC22))
# Autovalores e autovetores de b
EigB = eigen(PARESB)
# Coeficientes canônicos b
b = inv(sqrt(VC22)) * EigB.vectors
# Correlações canônicas
RUV = sqrt.(EigA.values)
# Ordenação em ordem decrescente
idx = sortperm(RUV, rev = true)
RUV = RUV[idx]
# Reordenar coeficientes a
a = a[:, idx]
# Primeira função canônica
a1 = a[:, 1]
b1 = b[:, sortperm(EigB.values, rev = true)[1]]
# Correlações entre variáveis canônicas e variáveis originais
DVC11 = diagm(diag(VC11))
DVC22 = diagm(diag(VC22))
IRDVC11 = inv(sqrt(DVC11))
IRDVC22 = inv(sqrt(DVC22))
RU1X1 = a1' * VC11 * IRDVC11 # U1 com X1
RU1X2 = a1' * VC12 * IRDVC22 # U1 com X2
RV1X1 = b1' * VC21 * IRDVC11 # V1 com X1
RV1X2 = b1' * VC22 * IRDVC22 # V1 com X2
# Dados transformados: variáveis canônicas
# As variáveis canônicas são combinações lineares das variáveis originais
# U = X1 * A (variáveis canônicas do primeiro conjunto)
# V = X2 * B (variáveis canônicas do segundo conjunto)
279
U = X1 * a # Matriz (n × p) de variáveis canônicas de X1
V = X2 * b # Matriz (n × q) de variáveis canônicas de X2
# Primeiro par de variáveis canônicas
U1 = U[:, 1] # Primeira variável canônica associada a X1
V1 = V[:, 1] # Primeira variável canônica associada a X2
# Verificação da correlação canônica
# Correlação entre o primeiro par canônico (U1, V1)
cor_U1V1 = cor(U1, V1)
# Estatística de Wilks e teste qui-quadrado
# V11* = V11 - V12 * V22^{-1} * V21
# V11 = VC11 (covariância de X1)
# V22 = VC22 (covariância de X2)
# V12 = VC12 (covariância cruzada)
# V21 = VC21 (transposta)
V11_star = VC11 - VC12 * inv(VC22) * VC21 
# Somas de quadrados 
SQTotal = det(VC11)
SQResiduo = det(V11_star)
# Critério de Wilks
L = SQResiduo / SQTotal
# Estatística qui-quadrado
p = size(X1, 2) # número de variáveis em X1
q = size(X2, 2) # número de variáveis em X2
chi2 = (-n + 0.5 * (p + q + 1)) * log(L)
df = (p - 1) * (q - 1) 
L # Estatística de Wilks
chi2 # Valor do qui-quadrado
_df = df # Graus de liberdade
p_valor_X = ccdf(Chisq(_df), chi2)
# Resumo dos resultados para as variáveis originais
# Conjunto X1 
# Linhas: variáveis originais de X1
# Colunas: correlação com U1 (RU1X1) e V1 (RV1X1)
Resumo_X1 = [RU1X1' RV1X1']
# Conjunto X2 
# Linhas: variáveis originais de X2
# Colunas: correlação com U1 (RU1X2) e V1 (RV1X2)
Resumo_X2 = [RU1X2' RV1X2']
280
# Análise de correlação a partir das variáveis padronizadas
# Padronização das variáveis
Z1 = (X1 .- mean(X1, dims=1)) ./ std(X1, dims=1)
Z2 = (X2 .- mean(X2, dims=1)) ./ std(X2, dims=1)
#Matrizes de correlação
R11 = cor(Z1) # Correlação de Z1
R22 = cor(Z2) # Correlação de Z2
R12 = cor(Z1, Z2) # Correlação cruzada
R21 = R12'
R = cor([Z1 Z2]) # Correlação total
# Determinação dos coeficientes canônicos padronizados
# Matriz associada aos coeficientes a*
PARESA_Z = inv(sqrt(R11)) * R12 * inv(R22) * R21 * inv(sqrt(R11))
EigA_Z = eigen(PARESA_Z)
a_Z = inv(sqrt(R11)) * EigA_Z.vectors
# Matriz associada aos coeficientes b*
PARESB_Z = inv(sqrt(R22)) * R21 * inv(R11) * R12 * inv(sqrt(R22))
EigB_Z = eigen(PARESB_Z)
b_Z = inv(sqrt(R22)) * EigB_Z.vectors
# Correlações canônicas padronizadas
RUV_Z = sqrt.(EigA_Z.values)
# Ordenação decrescente
idx_Z = sortperm(RUV_Z, rev=true)
RUV_Z = RUV_Z[idx_Z]
# Reordenar coeficientes
a_Z = a_Z[:, idx_Z]
# Primeira função canônica padronizada
a1_Z = a_Z[:, 1]
b1_Z = b_Z[:, sortperm(EigB_Z.values, rev=true)[1]]
# Correlações entre variáveis canônicas e variáveis originais padronizadas
RU1Z1 = a1_Z' * R11 # U1* com Z1
RU1Z2 = a1_Z' * R12 # U1* com Z2
RV1Z1 = b1_Z' * R21 # V1* com Z1
RV1Z2 = b1_Z' * R22 # V1* com Z2
# Variáveis canônicas padronizadas
U_Z = Z1 * a_Z
V_Z = Z2 * b_Z
U1_Z = U_Z[:, 1]
V1_Z = V_Z[:, 1]
# Verificação da correlação canônica padronizada
cor_U1Z_V1Z = cor(U1_Z, V1_Z)
281
# Estatística de Wilks (padronizada)
R11_star = R11 - R12 * inv(R22) * R21
SQTotal_Z = det(R11)
SQResiduo_Z = det(R11_star)
L_Z = SQResiduo_Z / SQTotal_Z
chi2_Z = (-n + 0.5 * (p + q + 1)) * log(L_Z)
df_Z = (p - 1) * (q - 1)
p_valor_Z = ccdf(Chisq(df_Z), chi2_Z)
# Resumo dos resultados — variáveis padronizadas
# Conjunto Z1
Resumo_Z1 = [
 RU1Z1' RV1Z1']
#Conjunto Z2
Resumo_Z2 = [RU1Z2' RV1Z2']