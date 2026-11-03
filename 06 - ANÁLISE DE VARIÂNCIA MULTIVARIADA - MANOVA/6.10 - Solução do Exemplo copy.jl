# SOLUÇÃO NO JULIA
#Matriz de dados geral (variáveis respostas)
X=[43.48 8.42 73.70 0.56 0.49
67.96 12.79 121.04 0.56 0.50
65.19 11.61 114.27 0.54 0.47
70.50 14.36 119.43 0.59 0.50
47.40 8.53 78.21 0.59 0.51
52.90 10.58 95.22 0.57 0.49
61.30 12.87 109.73 0.53 0.46
59.70 11.04 98.51 0.57 0.48
72.10 17.46 143.16 0.52 0.53
70.87 20.93 138.80 0.51 0.50
85.31 22.18 197.15 0.51 0.50
66.04 20.70 124.67 0.52 0.54
74.06 17.77 216.39 0.51 0.51
89.53 18.80 222.53 0.54 0.54
69.15 18.58 119.72 0.58 0.55
87.99 21.13 146.65 0.58 0.56
110.08 19.43 178.07 0.58 0.58
159.37 33.56 278.88 0.58 0.57
127.09 17.30 222.87 0.60 0.56
111.52 14.38 201.16 0.58 0.54
143.13 24.24 257.03 0.56 0.56
137.43 27.21 236.47 0.58 0.58
115.42 20.30 190.76 0.58 0.57
152.91 24.84 246.90 0.58 0.53]
#Matriz de dados para o grupo 1 (Tratamento 1)
X1=X[1:8,:]
#Matriz de dados para o grupo 2 (Tratamento 2)
X2=X[9:16,:]
#Matriz de dados para o grupo 3 (Tratamento 3)
X3=X[17:24,:]
#Tamanho da amostra do grupo 1 (n1) e número de variáveis (p)
(n1,p)=size(X1)
#Tamanho da amostra do grupo 2 (n2) e número de variáveis (p)
(n2,p)=size(X2)
#Tamanho da amostra do grupo 3 (n3) e número de variáveis (p)
(n3,p)=size(X3)
#Tamanho total das amostras
n=n1+n2+n3
#Vetor do tamanho das amostras dos grupos
ng=[n1;n2;n3]
192
#Número de grupos
g=size(ng,1)
#Vetor de médias para o grupo 1
MEDIAX1=mean(X1,dims=1)
#Vetor de médias para o grupo 2
MEDIAX2=mean(X2,dims=1)
#Vetor de médias para o grupo 3
MEDIAX3=mean(X3,dims=1)
#Vetor de médias dos grupos
MEDIAS=[MEDIAX1;MEDIAX2;MEDIAX3]
#Vetor de médias geral
MEDIAGERAL=mean(X,dims=1)
#Matriz de diferenças das médias dos grupos e a média geral
DIFMEDIA=MEDIAS.-MEDIAGERAL
#Matriz de diferenças das médias dos grupos e a média geral corrigida pelo tamanho dos grupos
DIFMEDIACOR=DIFMEDIA.*ng
#Matriz de somas de quadrados dos grupos
B=DIFMEDIACOR'*DIFMEDIA
#Matriz de diferenças dos valores de cada grupo e a média de seu respectivo 
grupo
DIFX1=X1.-MEDIAX1
DIFX2=X2.-MEDIAX2
DIFX3=X3.-MEDIAX3
DIF=[DIFX1;DIFX2;DIFX3]
#Matriz de somas de quadrados dos resíduos
W=DIF'*DIF
#Cálculo do valor de Lambda de Wilks
L=det(W)/det(B+W)
#Cálculo do valor de Lambda de Wilks corrigido
Lcorrigido=-(n-1-((p+g)/2))*log(L)
#Valor crítico de Qui-Quadrado 
QUIQUADRAINV= quantile(Chisq(p * (g - 1)), 0.99)
