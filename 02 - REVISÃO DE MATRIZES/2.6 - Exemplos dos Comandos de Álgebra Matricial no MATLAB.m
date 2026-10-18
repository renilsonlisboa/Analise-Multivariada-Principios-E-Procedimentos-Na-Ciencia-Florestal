%QUANTIDADE ESCALAR
A=2

%VETOR LINHA%
A=[2 3 4]
A =
 2 3 4
VETOR COLUNA
B=[2 ; 3; 4]
B =
 2
 3
 4
MATRIZ
A=[3 5 7 ; 8 9 10]
A =
 3 5 7
 8 9 10
ORDEM DE UMA MATRIZ
A=[5 6 7 ; 9 8 7 ;5 4 3]
A =
 5 6 7
 9 8 7
 5 4 3
[linhas , colunas]=size(A)
linhas =
 3
colunas =
 3
MATRIZ QUADRADA
A=[2 4 ; 5 6]
A =
40
 2 4
 5 6
MATRIZ NULIDADE OU ZERO
A=zeros(3,3)
A =
 0 0 0
 0 0 0
 0 0 0
DIAGONAL DE UMA MATRIZ
A=[5 6 7 ; 9 8 7 ; 5 4 3]
A =
 5 6 7
 9 8 7
 5 4 3
DIAGONAL=diag(A)
DIAGONAL =
 5
 8
 3
MATRIZ DIAGONAL
V=[3 ; 5 ; 1]
V =
 3
 5
 1
A=diag(V)
A =
 3 0 0
 0 5 0
 0 0 1
MATRIZ IDENTIDADE
A=eye(4)
A =
 1 0 0 0
 0 1 0 0
 0 0 1 0
 0 0 0 1
TRIÂNGULO SUPERIOR DA MATRIZ
A=[5 6 7 1 ; 9 8 7 3 ; 5 4 3 2 ; 7 9 5 1]
A =
 5 6 7 1
 9 8 7 3
 5 4 3 2
 7 9 5 1
B=tril(A)
B =
 5 6 7 1
 0 8 7 3
 0 0 3 2
 0 0 0 1
TRIÂNGULO INFERIOR DA MATRIZ
C=tril(A)
41
C =
 5 0 0 0
 9 8 0 0
 5 4 3 0
 7 9 5 1
MATRIZ SIMÉTRICA
A=[ 4 3 ; 3 5]
A =
 4 3
 3 5
MATRIZ ESCALAR
C=10*[1 2 ; 3 4]
C =
 10 20
 30 40
IGUALDADE DE MATRIZES
A=[4 2 ; 6 1]
A =
 4 2
 6 1
B=[4 2 ; 3 1]
B =
 4 2
 3 1
B=A
B =
 4 2
 6 1
TRANSPOSTA DE UMA MATRIZ
A=[4 5 6 ; 7 8 9]
A =
 4 5 6
 7 8 9
B=A'
B =
 4 7
 5 8
 6 9
ADIÇÃO DE MATRIZES
A=[2 2 ; 3 3]
A =
 2 2
 3 3
B=[ 5 6 ; 6 5]
B =
 5 6
 6 5
C=A+B
C =
 7 8
 9 8
42
SUBTRAÇÃO DE MATRIZES
D=A-B
D =
 -3 -4
 -3 -2
MULTIPLICAÇÃO DE VETORES
A=[2 5 -2 ]
A =
 2 5 -2
B=[3 ; -1 ; 4]
B =
 3
 -1
 4
C=A*B
C =
 -7
MULTIPLICAÇÃO DE MATRIZES
A=[2 2 ; 3 3]
A =
 2 2
 3 3
B=[5 6 ; 6 5]
B =
 5 6
 6 5
C = A * B
C = 
 22 22
 33 33
EQUAÇÕES LINEARES SIMULTÂNEAS NA FORMA MATRICIAL
A =[5 8; 9 4]
A =
 5 8
 9 4
B = [6;4]
B = 
 6
 4
X=inv(A)*B
X = 
 0.1538
 0.6538
B=A*X
B = 
 6.0000
 4.0000
DETERMINANTES
A=[ 3 2 ; 6 5]
A =
 3 2
 6 5
B=det(A)
43
B =
 3
SINGULARIDADE
A=[2 2 ; 2 2]
A =
 2 2
 2 2
B=det(A)
B =
 0
C=inv(A)
[Warning: Matrix is singular to working precision.] 
C =
 Inf Inf
 Inf Inf
MATRIZ INVERSA
A=[1 3 4 ; 5 6 7 ; 9 8 7]
A =
 1 3 4
 5 6 7
 9 8 7
B=inv(A)
B =
 -1.0000 0.7857 -0.2143
 2.0000 -2.0714 0.9286
 -1.0000 1.3571 -0.6429
TRAÇO DA MATRIZ
A=[1 3 4 ; 5 6 7 ; 9 8 7 ]
A =
 1 3 4
 5 6 7
 9 8 7
B=trace(A)
B =
 14
VALORES CARACTERÍSTICOS E VETORES CARACTERÍSITICOS (EIGENVALUES E 
EIGENVECTORS)
A=[ 1 2 ; 2 1]
A =
 1 2
 2 1
[e , l] = eig(A)
e =
 -0.7071 0.7071
 0.7071 0.7071
l =
 -1 0
 0 3
