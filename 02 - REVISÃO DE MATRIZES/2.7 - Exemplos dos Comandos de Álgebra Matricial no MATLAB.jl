#Pacotes
using LinearAlgebra
#Quantidade escalar
A=2
#Vetor linha
A=[2 3 4]
#Vetor coluna
B=[2;3;4]
B=A'
#Matriz
A=[3 5 7; 8 9 10]
#Ordem de uma matriz
A=[5 6 7; 9 8 7; 5 4 3]
size(A)
A_linhas, A_colunas = size(A)
A_linhas
A_colunas
#Matriz quadrada
A=[2 4; 5 6]
#Matriz nulidade ou zero
A = zeros(3,3)
#Diagonal de uma matriz
A=[5 6 7; 9 8 7; 5 4 3]
Diagonal=diag(A)
Diagonal=diagonal(A)
#Matriz diagonal
V = [3, 5, 1]
A=diagm(V)
A=Diagonal(V)
#Mariz identidade
A=Matrix(I, 4, 4)
#Triângulo superior da matriz
A=[5 6 7 1; 9 8 7 3; 5 4 3 2; 7 9 5 1]
B=UpperTriangular(A)
#Triângulo inferior da matriz
C=LowerTriangular(A)
45
#Matriz simétrica
A=[4 3; 1 5]
B=Symmetric(A)
#Matriz escalar
C=10*[1 2; 3 4]
#Igualdade de matrizes
A=[4 2; 6 1]
B=[4 2; 6 1]
B=A
#Transposta de uma matriz
A=[4 5 6; 7 8 9]
B=A'
#Adição de matrizes
A=[2 2; 3 3]
B=[5 6; 6 5]
C=A+B
#Subtração de matrizes
D=A-B
#Multiplicação de vetores
A=[2 5 -2]
B=[3;-3;4]
C=A*B
#Multiplicação de matrizes
A=[2 2; 3 3]
B=[5 6; 6 5]
C=A*B
#Equações lineares simultâneas na forma matricial
A=[5 8; 9 4]
B=[6;4]
X=inv(A)*B 
B=A*X
#Determinantes
A=[3 2; 6 5]
B=det(A)
#Singularidade
A=[2 2; 2 2]
B=det(A)
C=inv(A)
46
#Matriz adjunta
A=[4 5; 6 7]
B=Adjoint(A)
#Matriz inversa
A=[1 3 4; 5 6 7; 9 8 7]
B=inv(A)
#Traço da matriz
A=[1 3 4; 5 6 7; 9 8 7]
B=tr(A)
#Valores característicos e vetores característicos (eigenvalues e eigenvectors)
A=[1 2; 2 1]
B=eigen(A)