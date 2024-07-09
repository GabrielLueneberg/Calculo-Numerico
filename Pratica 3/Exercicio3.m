A = [15 2 6 8; 3 4 7 5; 7 7 8 9; 4 5 6 7];
n = size(A, 1);
Q = zeros(n, n);
R = zeros(n, n);

for j = 1 : n
    v = A(:, j);
    for i = 1 : j - 1
        R(i, j) = Q(:, i)' * A(:, j);
        v = v - R(i, j) * Q(:, i);
    end
    R(j, j) = norm(v);
    Q(:, j) = v / R(j, j);
end
disp('Matrizes Q e R após decomposição QR')
disp(" ")
Q
disp(" ")
R
disp(" ")
disp('Produto de Q e R')
disp(" ")
Q * R


%Autovalores pelo Metodo de Francis
disp('#############################################################################################################################################################################')
disp(" ")
disp('Autovalores pelo Método de Francis')
disp(" ")
tol = 1e-20;
erro = 1;
while (erro > tol)
    A = R * Q;
    [Q, R] = qr(A);
    erro = norm(tril(A, -1), 'fro');
end

for i = 1 : n
    autovalor = A(i, i)
end
disp(" ")
disp('Autovalores e autovetores pela função eig')
disp(" ")
[V,e] = eig(A)


%Autovalores pelo Metodo da potencia
disp('#############################################################################################################################################################################')
disp(" ")
disp('Autovalores pelo método da potência inversa')
disp(" ")
x = rand(n, 1);
lambda_ant = inf;
B = inv(A);
while true
    y = x / norm(x, 'fro');
    lambda = y' * B * y;
    x = B * y;
    if abs(lambda - lambda_ant) < tol
        break;
    end
    lambda_ant = lambda;
end
lambda = 1/lambda
disp(" ")
y
disp(" ")
disp('Autovalores e autovetores pela função eig')
disp(" ")
[V,e] = eig(A)
