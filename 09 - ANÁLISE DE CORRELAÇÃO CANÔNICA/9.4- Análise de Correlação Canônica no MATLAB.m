% ANÁLISE DE CORRELAÇÃO CANÔNICA
% MATRIZ DE DADOS X (N × (P+Q))
X = [
    174.45 31.64 206.09 19.74 20.2 7.7 0.65
    84.41 8.76 93.17 15.12 19.1 4.9 0.5
    60.97 4.91 65.88 12.89 17.3 3 0.57
    47.25 2.06 49.3 11.62 16.9 7.2 0.52
    114.06 27.95 142 16.55 20.5 11.8 0.56
    157.39 37.48 194.87 19.1 21.9 7 0.6
    72.03 15.14 87.17 14.32 18.4 9.2 0.56
    88.98 15.66 104.63 15.12 19.3 12.1 0.58
    45.29 5.96 51.25 11.94 16.7 9.5 0.5
    101.2 18.46 119.65 15.28 19.5 7.7 0.59
    94.16 12.22 106.38 16.23 18.3 4.6 0.52
    83.79 10.13 93.92 15.12 19.3 3.8 0.59
    200 54.52 254.51 20.53 19.5 9.5 0.58
    49.66 6.04 55.7 13.21 18.4 2.8 0.57
    20.96 3.59 24.55 8.75 14.4 11.2 0.48
    71.9 15.34 87.24 13.05 18 3.5 0.62
    115.57 21.25 136.82 16.23 19.6 9 0.62
    69.3 11.5 80.8 12.25 17.9 3.8 0.63
    16.46 2.49 18.94 7.64 9.9 6.7 0.51
    107.51 19.67 127.17 15.44 19.2 6.8 0.56
];

% DEFINIÇÃO DOS CONJUNTOS DE VARIÁVEIS
% X1: PRIMEIRO CONJUNTO (N × P)
X1 = X(:,1:3)
% X2: SEGUNDO CONJUNTO (n × q)
X2 = X(:,4:7)
% NÚMERO DE OBSERVAÇÕES
n = size(X1,1)
% MATRIZES DE VARIÂNCIA–COVARIÂNCIA
% Calcula a matriz de covariância conjunta de todas as variáveis
% Dimensão: (p+q) × (p+q)
% As p primeiras colunas/linhas correspondem a X1
% As q últimas colunas/linhas correspondem a X2
C = cov([X1 X2])
% NÚMERO DE VARIÁVEIS EM CADA CONJUNTO
p = size(X1,2); % quantidade de variáveis no primeiro conjunto
q = size(X2,2); % quantidade de variáveis no segundo conjunto
% MATRIZ DE COVARIÂNCIA DO PRIMEIRO CONJUNTO DE VARIÁVEIS
% Dimensão: p × p
VC11 = C(1:p, 1:p)
% MATRIZ DE COVARIÂNCIA DO SEGUNDO CONJUNTO DE VARIÁVEIS
% Dimensão: q × q
VC22 = C(p+1:p+q, p+1:p+q)
% MATRIZ DE COVARIÂNCIA CRUZADA ENTRE X1 E X2
% Cada elemento representa cov(X1_i, X2_j)
% Dimensão: p × q
VC12 = C(1:p, p+1:p+q)
% TRANSPOSTA DA MATRIZ DE COVARIÂNCIA CRUZADA
% Representa cov(X2, X1)
% Dimensão: q × p
VC21 = VC12'
% MATRIZ DE COVARIÂNCIA TOTAL
% Dimensão: (p+q) × (p+q)
VC = C
% DETERMINAÇÃO DOS COEFICIENTES CANÔNICOS
% MATRIZ ASSOCIADA AOS COEFICIENTES A
PARESA = inv(sqrtm(VC11)) * VC12 * inv(VC22) * VC21 * inv(sqrtm(VC11))
% AUTOVALORES E AUTOVETORES DE A
[VA,DA] = eig(PARESA)
% COEFICIENTES CANÔNICOS A
a = inv(sqrtm(VC11)) * VA
% MATRIZ ASSOCIADA AOS COEFICIENTES B
PARESB = inv(sqrtm(VC22)) * VC21 * inv(VC11) * VC12 * inv(sqrtm(VC22))
% AUTOVALORES E AUTOVETORES DE B
[VB,DB] = eig(PARESB)
% COEFICIENTES CANÔNICOS B
b = inv(sqrtm(VC22)) * VB
% CORRELAÇÕES CANÔNICAS
% As correlações canônicas são as raízes dos autovalores
RUV = sqrt(diag(DA))
% ORDENAÇÃO EM ORDEM DECRESCENTE
[RUV,idx] = sort(RUV,'descend')
a = a(:,idx)
% PRIMEIRA FUNÇÃO CANÔNICA
a1 = a(:,1)
[~,idxB] = sort(diag(DB),'descend');
b1 = b(:,idxB(1));
% CORRELAÇÕES ENTRE VARIÁVEIS CANÔNICAS E VARIÁVEIS ORIGINAIS
DVC11 = diag(diag(VC11))
DVC22 = diag(diag(VC22))
IRDVC11 = inv(sqrtm(DVC11))
IRDVC22 = inv(sqrtm(DVC22))
RU1X1 = a1' * VC11 * IRDVC11 % U1 com X1
RU1X2 = a1' * VC12 * IRDVC22 % U1 com X2
RV1X1 = b1' * VC21 * IRDVC11 % V1 com X1
RV1X2 = b1' * VC22 * IRDVC22 % V1 com X2
% DADOS TRANSFORMADOS: VARIÁVEIS CANÔNICAS
U = X1 * a % Variáveis canônicas do primeiro conjunto
V = X2 * b % Variáveis canônicas do segundo conjunto
U1 = U(:,1) % Primeira variável canônica de X1
V1 = V(:,1) % Primeira variável canônica de X2
% VERIFICAÇÃO DA CORRELAÇÃO CANÔNICA
cor_U1V1 = corr(U1,V1)
% ESTATÍSTICA DE WILKS E TESTE QUI-QUADRADO
% MATRIZ DE RESÍDUOS DA REGRESSÃO
V11_star = VC11 - VC12 * inv(VC22) * VC21
% SOMAS DE QUADRADOS (DETERMINANTES)
SQTotal = det(VC11)
SQResiduo = det(V11_star)
% CRITÉRIO DE WILKS
L = SQResiduo / SQTotal
% ESTATÍSTICA QUI-QUADRADO
p = size(X1,2)
q = size(X2,2)
chi2 = (-n + 0.5*(p+q+1)) * log(L)
df = (p-1)*(q-1)
p_valor_X = 1 - chi2cdf(chi2,df)
% RESUMO DOS RESULTADOS — VARIÁVEIS ORIGINAIS
Resumo_X1 = [RU1X1' RV1X1']
Resumo_X2 = [RU1X2' RV1X2']
% ANÁLISE DE CORRELAÇÃO A PARTIR DAS VARIÁVEIS PADRONIZADAS
% Padronização das variáveis
Z1 = (X1 - mean(X1)) ./ std(X1);
Z2 = (X2 - mean(X2)) ./ std(X2);
% MATRIZES DE CORRELAÇÃO
R11 = corr(Z1)
R22 = corr(Z2)
R12 = corr(Z1,Z2)
R21 = R12'
R = corr([Z1 Z2])
% COEFICIENTES CANÔNICOS PADRONIZADOS
PARESA_Z = inv(sqrtm(R11)) * R12 * inv(R22) * R21 * inv(sqrtm(R11))
[VAZ,DAZ] = eig(PARESA_Z)
a_Z = inv(sqrtm(R11)) * VAZ
PARESB_Z = inv(sqrtm(R22)) * R21 * inv(R11) * R12 * inv(sqrtm(R22))
[VBZ,DBZ] = eig(PARESB_Z)
b_Z = inv(sqrtm(R22)) * VBZ
% CORRELAÇÕES CANÔNICAS PADRONIZADAS
RUV_Z = sqrt(diag(DAZ))
[RUV_Z,idx_Z] = sort(RUV_Z,'descend')
a_Z = a_Z(:,idx_Z)
% PRIMEIRA FUNÇÃO CANÔNICA PADRONIZADA
a1_Z = a_Z(:,1)
[~,idxBZ] = sort(diag(DBZ),'descend');
b1_Z = b_Z(:,idxBZ(1))
% CORRELAÇÕES ENTRE VARIÁVEIS CANÔNICAS E VARIÁVEIS PADRONIZADAS
RU1Z1 = a1_Z' * R11
RU1Z2 = a1_Z' * R12
RV1Z1 = b1_Z' * R21
RV1Z2 = b1_Z' * R22
% VARIÁVEIS CANÔNICAS PADRONIZADAS
U_Z = Z1 * a_Z
V_Z = Z2 * b_Z
U1_Z = U_Z(:,1)
V1_Z = V_Z(:,1)
% VERIFICAÇÃO DA CORRELAÇÃO CANÔNICA PADRONIZADA
cor_U1Z_V1Z = corr(U1_Z,V1_Z)
% ESTATÍSTICA DE WILKS (PADRONIZADA)
R11_star = R11 - R12 * inv(R22) * R21
SQTotal_Z = det(R11)
SQResiduo_Z = det(R11_star)
L_Z = SQResiduo_Z / SQTotal_Z
chi2_Z = (-n + 0.5*(p+q+1)) * log(L_Z)
df_Z = (p-1)*(q-1)
p_valor_Z = 1 - chi2cdf(chi2_Z,df_Z)
% RESUMO — VARIÁVEIS PADRONIZADAS
Resumo_Z1 = [RU1Z1' RV1Z1']
Resumo_Z2 = [RU1Z2' RV1Z2']
