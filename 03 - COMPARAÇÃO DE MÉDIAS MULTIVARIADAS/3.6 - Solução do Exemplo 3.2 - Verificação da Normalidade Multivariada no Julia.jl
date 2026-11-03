#SOLUÇÃO DO EXEMPLO 3.2 NO JULIA
#Ordenação das distâncias generalizadas cáculadas no exemplo 3.1
D2JORD=sort(D2j,dims=1)
#Obtenção do vetor de posições (j) 
POS=collect(1:1:n)
#Valores de alfa para cada posição j
ALFA=(POS.-0.5)./n
#Valores de Qui-quadrado para cada posição j
y = quantile(Chisq(p), ALFA)
#Gráfico para avaliação da normalidade p variada
plot(D2JORD,y)
gcf()
#Ajuste da equação da reta 
lr = linregress(D2JORD, y)
#Estimativa dos valores da reta para o gráfico
ybest = lr(D2JORD)
#Gráfico para avaliação da normalidade p variada
plot(D2JORD,ybest)
xlabel("DISTÂNCIA GENERALIZADA")
ylabel("VALOR DE QUI-QUADRADO")
title("AVALIAÇÃO DA NORMALIDADE p VARIADA")
gcf()