A = [15 2 6 8; 3 4 7 5; 7 7 8 9; 4 5 6 7]
I = eye(size(A));
n = size(A, 1);
[Q, R] = qr(A);

tol = 1e-20;
erro = 1;
while (erro > tol)
    A = R * Q;
    [Q, R] = qr(A);
    erro = norm(tril(A, -1), 'fro');
end


for i = 1 : n
    autovalor = A(i, i)
    autovetor = null(A - autovalor * I)
    disp(" ")

end



disp(" ")
disp('Autovalores e autovetores pela função eig')
disp(" ")
[V,e] = eig(A)











