%% ANÁLISE MULTIVARIADA - TESTE DE LAMBDA DE WILKS
% SOLUÇÃO NO MATLAB

%% 1. MATRIZ GERAL DE DADOS (VARIÁVEIS RESPOSTAS)

X = [
    43.4800  8.4200  73.7000  0.5600  0.4900
    67.9600 12.7900 121.0400  0.5600  0.5000
    65.1900 11.6100 114.2700  0.5400  0.4700
    70.5000 14.3600 119.4300  0.5900  0.5000
    47.4000  8.5300  78.2100  0.5900  0.5100
    52.9000 10.5800  95.2200  0.5700  0.4900
    61.3000 12.8700 109.7300  0.5300  0.4600
    59.7000 11.0400  98.5100  0.5700  0.4800
    72.1000 17.4600 143.1600  0.5200  0.5300
    70.8700 20.9300 138.8000  0.5100  0.5000
    85.3100 22.1800 197.1500  0.5100  0.5000
    66.0400 20.7000 124.6700  0.5200  0.5400
    74.0600 17.7700 216.3900  0.5100  0.5100
    89.5300 18.8000 222.5300  0.5400  0.5400
    69.1500 18.5800 119.7200  0.5800  0.5500
    87.9900 21.1300 146.6500  0.5800  0.5600
    110.0800 19.4300 178.0700 0.5800  0.5800
    159.3700 33.5600 278.8800 0.5800  0.5700
    127.0900 17.3000 222.8700 0.6000  0.5600
    111.5200 14.3800 201.1600 0.5800  0.5400
    143.1300 24.2400 257.0300 0.5600  0.5600
    137.4300 27.2100 236.4700 0.5800  0.5800
    115.4200 20.3000 190.7600 0.5800  0.5700
    152.9100 24.8400 246.9000 0.5800  0.5300
];

%% 2. MATRIZES DE DADOS POR GRUPO (TRATAMENTOS)

% Grupo 1 (Tratamento 1)
X1 = X(1:8,:);

% Grupo 2 (Tratamento 2)
X2 = X(9:16,:);

% Grupo 3 (Tratamento 3)
X3 = X(17:24,:);

%% 3. TAMANHO DAS AMOSTRAS E NÚMERO DE VARIÁVEIS

[n1,p] = size(X1);
[n2,p] = size(X2);
[n3,p] = size(X3);

% Tamanho total das amostras
n = n1 + n2 + n3;

% Vetor com o tamanho das amostras dos grupos
ng = [n1; n2; n3];

% Número de grupos
g = size(ng,1);

%% 4. VETORES DE MÉDIAS DOS GRUPOS

MEDIAX1 = mean(X1);
MEDIAX2 = mean(X2);
MEDIAX3 = mean(X3);

% Matriz de médias dos grupos
MEDIAS = [
    MEDIAX1
    MEDIAX2
    MEDIAX3
];

% Vetor de médias gerais
MEDIAGERAL = mean(X);

%% 5. MATRIZ DE SOMA DE QUADRADOS ENTRE GRUPOS

% Diferenças entre as médias dos grupos e a média geral
DIFMEDIA = MEDIAS - MEDIAGERAL;

% Correção pelo tamanho das amostras dos grupos
DIFMEDIACOR = DIFMEDIA .* ng;

% Matriz de soma de quadrados entre grupos
B = DIFMEDIACOR' * DIFMEDIA;

%% 6. MATRIZ DE SOMA DE QUADRADOS DOS RESÍDUOS

% Diferenças entre os valores observados e as médias dos respectivos grupos
DIFX1 = X1 - MEDIAX1;
DIFX2 = X2 - MEDIAX2;
DIFX3 = X3 - MEDIAX3;

% Matriz geral de diferenças
DIF = [
    DIFX1
    DIFX2
    DIFX3
];

% Matriz de soma de quadrados dos resíduos
W = DIF' * DIF;

%% 7. CÁLCULO DO LAMBDA DE WILKS

L = det(W) / det(B + W);

%% 8. CÁLCULO DO LAMBDA DE WILKS CORRIGIDO

LCORRIGIDO = -(n - 1 - ((p + g)/2)) * log(L);

%% 9. VALOR CRÍTICO DA DISTRIBUIÇÃO QUI-QUADRADO

QUI2 = chi2inv(0.99, p * (g - 1));