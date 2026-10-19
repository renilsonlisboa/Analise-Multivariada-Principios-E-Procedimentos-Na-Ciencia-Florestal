% ANÁLISE DE AGRUPAMENTO COM DISTÂNCIA EUCLIDIANA E LIGAÇÃO SIMPLES

% DADOS DE ENTRADA
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
 1.1800 0.3200 -0.3200 -0.7700
 -0.1700 1.7600 -1.5600 -0.9600
 0.9700 0.8000 -0.8100 -0.9600
];

% NÚMERO DE CASOS
numero_casos = size(X,1);

% MATRIZ DE DISTÂNCIAS EUCLIDIANAS
matriz_distancias = zeros(numero_casos);

for i = 1:numero_casos
    for j = i+1:numero_casos
        matriz_distancias(i,j) = sqrt(sum((X(i,:) - X(j,:)).^2));
        matriz_distancias(j,i) = matriz_distancias(i,j);
    end
end

% INICIALIZAÇÃO DOS CLUSTERS
conglomerados = num2cell(1:numero_casos);
distancias_atuais = matriz_distancias;
matriz_ligacao = [];
matriz_cofenetica = zeros(numero_casos);

% IDENTIFICAÇÃO DOS CONGLOMERADOS
identificador_conglomerado = numero_casos;

% ALGORITMO HIERÁRQUICO – LIGAÇÃO SIMPLES
while numel(conglomerados) > 1

    matriz_temporaria = distancias_atuais + ...
        diag(inf(size(distancias_atuais,1),1));

    [distancia_minima, indice] = min(matriz_temporaria(:));
    [i,j] = ind2sub(size(matriz_temporaria), indice);

    conglomerado_A = conglomerados{i};
    conglomerado_B = conglomerados{j};

    % ATUALIZAÇÃO DA MATRIZ COFENÉTICA
    for p = conglomerado_A
        for q = conglomerado_B
            if matriz_cofenetica(p,q) == 0
                matriz_cofenetica(p,q) = distancia_minima;
                matriz_cofenetica(q,p) = distancia_minima;
            end
        end
    end

    % ATUALIZAÇÃO DA MATRIZ DE LIGAÇÃO
    identificador_conglomerado = identificador_conglomerado + 1;

    matriz_ligacao = [matriz_ligacao; ...
        conglomerado_A(1), conglomerado_B(1), distancia_minima];

    novo_conglomerado = [conglomerado_A, conglomerado_B];

    % REMOÇÃO DOS CONGLOMERADOS UNIFICADOS
    indices_remocao = sort([i j], 'descend');
    conglomerados(indices_remocao) = [];
    distancias_atuais(indices_remocao,:) = [];
    distancias_atuais(:,indices_remocao) = [];

    % DISTÂNCIAS DO NOVO CONGLOMERADO
    numero_conglomerados = numel(conglomerados);
    novas_distancias = zeros(1,numero_conglomerados);

    for k = 1:numero_conglomerados
        lista_distancias = [];

        for a = novo_conglomerado
            for b = conglomerados{k}
                lista_distancias(end+1) = matriz_distancias(a,b);
            end
        end

        novas_distancias(k) = min(lista_distancias);
    end

    distancias_atuais = [distancias_atuais, novas_distancias';
                         novas_distancias, 0];

    conglomerados{end+1} = novo_conglomerado;
end

% CÁLCULO DA CORRELAÇÃO COFENÉTICA

% VETOR DAS DISTÂNCIAS ORIGINAIS
distancias_originais = ...
    matriz_distancias(triu(true(numero_casos),1));

% VETOR DAS DISTÂNCIAS COFENÉTICAS
distancias_cofeneticas = ...
    matriz_cofenetica(triu(true(numero_casos),1));

% CORRELAÇÃO COFENÉTICA DE PEARSON
correlacao_cofenetica = ...
    corr(distancias_originais, distancias_cofeneticas);

% RESULTADOS
disp('Correlação cofenética:');
disp(correlacao_cofenetica);

% DENDROGRAMA
figure;
dendrogram(matriz_ligacao);
title('Dendrograma - Ligação Simples');

% MAIOR DISTÂNCIA NO DENDROGRAMA
maior_distancia = max(matriz_ligacao(:,3));

% FUNÇÕES PARA TESTE DE DISTÂNCIAS E LIGAÇÕES

% TIPOS DE DISTÂNCIA
dist_types = { ...
    'euclidean', ...
    'seuclidean', ...
    'cityblock', ...
    'mahalanobis', ...
    'minkowski'};

% TIPOS DE LIGAÇÃO
link_types = { ...
    'single', ...    % vizinho mais próximo
    'complete', ...  % vizinho mais distante
    'average', ...   % média
    'centroid', ...  % centróide
    'ward'};         % Ward

% MATRIZ PARA ARMAZENAR RESULTADOS
R = zeros(length(dist_types), length(link_types));

% LOOP PRINCIPAL PARA COMBINAR DISTÂNCIAS E MÉTODOS
for i = 1:length(dist_types)

    % VETOR DE DISTÂNCIAS
    Y = pdist(X, dist_types{i});

    for j = 1:length(link_types)

        % WARD SÓ É MATEMATICAMENTE CONSISTENTE COM DISTÂNCIA EUCLIDIANA
        if strcmp(link_types{j}, 'ward') && ...
                ~strcmp(dist_types{i}, 'euclidean') && ...
                ~strcmp(dist_types{i}, 'seuclidean')

            R(i,j) = NaN;
            continue
        end

        % AGRUPAMENTO HIERÁRQUICO
        Z = linkage(Y, link_types{j});

        % CORRELAÇÃO COFENÉTICA
        R(i,j) = cophenet(Z, Y);
    end
end

% RESULTADOS DA CORRELAÇÃO COFENÉTICA PARA COMBINAÇÕES
% DE TIPOS DE DISTÂNCIAS E MÉTODOS DE LIGAÇÃO
T = array2table(R, ...
    'RowNames', dist_types, ...
    'VariableNames', link_types);

disp(T);

% MELHOR COMBINAÇÃO – DISTÂNCIA EUCLIDIANA AO QUADRADO,
% COM LIGAÇÃO CENTRÓIDE

% DENDROGRAMA
Z = linkage(Y, 'centroid');
dendrogram(Z, 0);

% RESULTADOS PARA VARIÁVEIS ORIGINAIS (NÃO-PADRONIZADAS)
disp(T);