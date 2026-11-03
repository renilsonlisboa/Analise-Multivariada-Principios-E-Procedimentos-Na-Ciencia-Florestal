%% ANÁLISE DISCRIMINANTE MULTIVARIADA
% EXEMPLO 7.2 E CLASSIFICAÇÃO (EXEMPLO 7.3)

%% 1. MATRIZES DE DADOS DOS GRUPOS

% Grupo 1
X1 = [
    0.47 68.62 28.10 5433
    0.47 69.02 27.85 5235
    0.48 69.10 27.67 5092.2
    0.51 69.10 27.67 5092.2
    0.53 69.13 28.23 5312
    0.55 69.10 27.67 5092.2
    0.45 69.10 27.67 5092.2
    0.50 69.10 27.67 5092.2
    0.54 71.19 25.52 4991
    0.57 70.68 27.19 5018
];

% Grupo 2
X2 = [
    0.49 73.44 24.92 4694
    0.51 74.27 23.93 4759
    0.49 74.46 23.53 4870
    0.57 72.28 24.97 4691
    0.56 72.37 25.12 4776
    0.55 70.23 26.77 4878
    0.52 73.49 23.83 4724
    0.56 72.39 24.89 4815
    0.54 73.16 24.37 4705
    0.51 72.62 23.95 4855
];

% Grupo 3
X3 = [
    0.56 80.99 18.50 4434
    0.58 74.21 23.69 4442
    0.51 73.77 24.34 4705
    0.64 76.84 21.61 4613
    0.54 74.01 23.99 4682
    0.52 72.19 25.79 4708
    0.34 71.55 25.49 4683
    0.59 73.32 24.41 4653.8
    0.52 77.53 21.59 4607
    0.58 74.73 23.30 4607
];

%% 2. MATRIZ GERAL DE DADOS

X = [X1; X2; X3];

% Número de grupos
g = 3;

%% 3. PADRONIZAÇÃO DOS DADOS

% Média e desvio-padrão das variáveis
mediaX = mean(X);
desvX = std(X);

% Dados padronizados de cada grupo
Z1 = (X1 - mediaX) ./ desvX;
Z2 = (X2 - mediaX) ./ desvX;
Z3 = (X3 - mediaX) ./ desvX;

%% 4. VETORES MÉDIOS DOS GRUPOS

MIX1 = mean(Z1);
MIX2 = mean(Z2);
MIX3 = mean(Z3);

%% 5. TAMANHO DAS AMOSTRAS E NÚMERO DE VARIÁVEIS

[n1,p] = size(X1);
[n2,p] = size(X2);
[n3,p] = size(X3);

% Vetor médio geral
MI = (n1*MIX1 + n2*MIX2 + n3*MIX3) / (n1+n2+n3);

%% 6. MATRIZES DE VARIÂNCIA-COVARIÂNCIA

S1 = cov(Z1);
S2 = cov(Z2);
S3 = cov(Z3);

%% 7. MATRIZ DE VARIÂNCIA-COVARIÂNCIA COMBINADA

Sp = ((n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3) / ...
     (n1+n2+n3-g);

%% 8. MATRIZ DE VARIAÇÃO ENTRE GRUPOS (B)

B = n1*(MIX1-MI)'*(MIX1-MI) + ...
    n2*(MIX2-MI)'*(MIX2-MI) + ...
    n3*(MIX3-MI)'*(MIX3-MI);

%% 9. MATRIZ DE VARIAÇÃO DENTRO DOS GRUPOS (W)

W = (n1+n2+n3-g) * Sp;

%% 10. AUTOVALORES DE W^(-1)B

WINV = inv(W);

WINVB = WINV * B;

%% 11. NÚMERO DE FUNÇÕES DISCRIMINANTES

s = min(g-1,p);

%% 12. AUTOVALORES E AUTOVETORES

[Vec,Val] = eig(WINVB);

% Extração dos autovalores
l = diag(Val);

% Ordenação dos autovalores em ordem decrescente
[l,ordem] = sort(l,'descend');

% Ordenação correspondente dos autovetores
Vec = Vec(:,ordem);

% Seleção dos autovalores correspondentes às funções discriminantes
l = l(1:s);

%% 13. AUTOVALORES ACUMULADOS

AVA = cumsum(l);

%% 14. VARIÂNCIA EXPLICADA ACUMULADA

VE = AVA / AVA(end) * 100;

%% 15. CORRELAÇÃO CANÔNICA

r = sqrt(abs(l) ./ (1 + abs(l)));

%% 16. TESTE LAMBDA DE WILKS

% Lambda de Wilks para a primeira função discriminante
Wilks1 = prod(1 ./ (1+l));

% Lambda de Wilks para as funções discriminantes restantes
Wilks = 1 ./ (1 + l(2:end));

%% 17. TESTE QUI-QUADRADO PARA AS FUNÇÕES DISCRIMINANTES

% Número de grupos
q = g;

% Número total de observações
n = n1+n2+n3;

% Valores de qui-quadrado
Qui2W1 = -(n - 1 - 0.5*(p+q)) * log(Wilks1);

Qui2W2 = -(n - 1 - 0.5*(p+q)) * log(Wilks);

%% 18. GRAUS DE LIBERDADE

% Primeira função discriminante
m = 1;
gl1 = (p-m+1)*(q-m);

% Segunda função discriminante
m = 2;
gl2 = (p-m+1)*(q-m);

%% 19. NÍVEIS DESCRITIVOS DOS TESTES

Alpha1 = 1 - chi2cdf(Qui2W1,gl1);

Alpha2 = 1 - chi2cdf(Qui2W2,gl2);

%% 20. OBTENÇÃO DOS COEFICIENTES DAS FUNÇÕES DISCRIMINANTES

% Cálculo de B0 (sem n)
B0 = (MIX1-MI)'*(MIX1-MI) + ...
     (MIX2-MI)'*(MIX2-MI) + ...
     (MIX3-MI)'*(MIX3-MI);

% Matriz inversa de W multiplicada por B0
WINVB0 = WINV * B0;

% Autovalores e autovetores
[eigVec,eigVal] = eig(WINVB0);

% Ordenação dos autovalores em ordem decrescente
[l,ordem] = sort(diag(eigVal),'descend');

% Ordenação correspondente dos autovetores
eigVec = eigVec(:,ordem);

% Seleção dos autovetores correspondentes às funções
e = eigVec(:,1:s);

%% 21. FATORES DE ESCALONAMENTO E AUTOVETORES ESCALONADOS

% Primeiro autovetor
e1 = e(:,1);

% Segundo autovetor
e2 = e(:,2);

% Fatores de escalonamento
k1 = 1/sqrt(e1' * Sp * e1);
k2 = 1/sqrt(e2' * Sp * e2);

% Autovetores escalonados
e1escal = k1 * e1;
e2escal = k2 * e2;

% Matriz dos coeficientes não padronizados
ee = [e1escal e2escal];

%% 22. COEFICIENTES DAS FUNÇÕES DISCRIMINANTES PADRONIZADAS

D = diag(Sp);

DR = sqrt(D);

eep = DR .* ee;

%% 23. MATRIZ DE ESTRUTURAS

Spd = diag(DR);

Rp = inv(Spd) * Sp * inv(Spd);

%% 24. MATRIZES DE ESTRUTURA DAS FUNÇÕES DISCRIMINANTES

% Primeira função discriminante
c1 = e1escal' * Spd;

me1 = Rp * c1';

% Segunda função discriminante
c2 = e2escal' * Spd;

me2 = Rp * c2';

% Matriz de estrutura final
me = [me1 me2];

%% 25. ESCORES DAS FUNÇÕES DISCRIMINANTES PADRONIZADAS

% Grupo 1
G1y1 = Z1 * eep(:,1);
G1y2 = Z1 * eep(:,2);

% Grupo 2
G2y1 = Z2 * eep(:,1);
G2y2 = Z2 * eep(:,2);

% Grupo 3
G3y1 = Z3 * eep(:,1);
G3y2 = Z3 * eep(:,2);

%% 26. GRÁFICO DAS FUNÇÕES DISCRIMINANTES

figure;
hold on;

scatter(G1y1,G1y2,'filled');
scatter(G2y1,G2y2,'filled');
scatter(G3y1,G3y2,'filled');

legend('Grupo 1','Grupo 2','Grupo 3');
xlabel('Função discriminante 1');
ylabel('Função discriminante 2');

grid on;
hold off;

%% 27. CENTRÓIDES DE CADA GRUPO

CG1 = [mean(G1y1) mean(G1y2)];
CG2 = [mean(G2y1) mean(G2y2)];
CG3 = [mean(G3y1) mean(G3y2)];

%% 28. CLASSIFICAÇÃO DAS OBSERVAÇÕES - EXEMPLO 7.3

% DISTÂNCIAS DO GRUPO 1 AOS CENTRÓIDES

dG1C1 = sqrt((G1y1-CG1(1)).^2 + (G1y2-CG1(2)).^2);
dG1C2 = sqrt((G1y1-CG2(1)).^2 + (G1y2-CG2(2)).^2);
dG1C3 = sqrt((G1y1-CG3(1)).^2 + (G1y2-CG3(2)).^2);

% Matriz de distâncias do grupo 1
DG1 = [dG1C1 dG1C2 dG1C3];

% Classificação pelo centróide mais próximo
[DG1,GrupoClassificado1] = min(DG1,[],2);

% Contagem de observações classificadas em cada grupo
contagem1 = histcounts(GrupoClassificado1,1:4);

%% 29. DISTÂNCIAS E CLASSIFICAÇÃO DO GRUPO 2

dG2C1 = sqrt((G2y1-CG1(1)).^2 + (G2y2-CG1(2)).^2);
dG2C2 = sqrt((G2y1-CG2(1)).^2 + (G2y2-CG2(2)).^2);
dG2C3 = sqrt((G2y1-CG3(1)).^2 + (G2y2-CG3(2)).^2);

% Matriz de distâncias do grupo 2
DG2 = [dG2C1 dG2C2 dG2C3];

% Classificação pelo centróide mais próximo
[DG2,GrupoClassificado2] = min(DG2,[],2);

% Contagem de observações classificadas em cada grupo
contagem2 = histcounts(GrupoClassificado2,1:4);

%% 30. DISTÂNCIAS E CLASSIFICAÇÃO DO GRUPO 3

dG3C1 = sqrt((G3y1-CG1(1)).^2 + (G3y2-CG1(2)).^2);
dG3C2 = sqrt((G3y1-CG2(1)).^2 + (G3y2-CG2(2)).^2);
dG3C3 = sqrt((G3y1-CG3(1)).^2 + (G3y2-CG3(2)).^2);

% Matriz de distâncias do grupo 3
DG3 = [dG3C1 dG3C2 dG3C3];

% Classificação pelo centróide mais próximo
[DG3,GrupoClassificado3] = min(DG3,[],2);

% Contagem de observações classificadas em cada grupo
contagem3 = histcounts(GrupoClassificado3,1:4);

%% 31. MATRIZ DE CONFUSÃO

% Grupos de origem
GrupoReal = [
    ones(n1,1)
    2*ones(n2,1)
    3*ones(n3,1)
];

% Grupos classificados
GrupoPredito = [
    GrupoClassificado1
    GrupoClassificado2
    GrupoClassificado3
];

% Matriz de confusão
MatrizConfusao = confusionmat(GrupoReal,GrupoPredito);

%% 32. CLASSIFICAÇÃO DE UMA NOVA OBSERVAÇÃO

Xnovo = [0.62 76.78 21.62 4614.11];

% Padronização da nova observação
Znovo = (Xnovo - mediaX) ./ desvX;

% Escores das funções discriminantes padronizadas
Ny1 = Znovo * eep(:,1);
Ny2 = Znovo * eep(:,2);

%% 33. DISTÂNCIAS DA NOVA OBSERVAÇÃO AOS CENTRÓIDES

dNovoC1 = sqrt((Ny1-CG1(1)).^2 + (Ny2-CG1(2)).^2);
dNovoC2 = sqrt((Ny1-CG2(1)).^2 + (Ny2-CG2(2)).^2);
dNovoC3 = sqrt((Ny1-CG3(1)).^2 + (Ny2-CG3(2)).^2);

% Matriz de distâncias
DNovo = [dNovoC1 dNovoC2 dNovoC3];

% Classificação da nova observação
[DNovo,GrupoClassificadonovo] = min(DNovo,[],2);