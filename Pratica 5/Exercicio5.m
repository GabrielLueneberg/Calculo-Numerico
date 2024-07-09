p = 4;
D = load("-ascii", "Table31.txt");
[m, n] = size(D);

y = D(:, 1);

A = zeros(m, p+1);

for j = 1:m
  for i = 1:p+1
    A(j, i) = y(j)^(i-1);
  end
end

disp('Matrix A:');
disp(A);

b = D(:, 2)

xa = (A' * A) \ (A' * b)


xb = pinv(A) * b

valores_eixo_x = linspace(min(y), max(y), 100);
valores_eixo_y = zeros(size(valores_eixo_x));

for i = 1:p+1
  valores_eixo_y += xb(i) * (valores_eixo_x).^(i-1);
end

plot(y, b, "r.", 'MarkerSize', 10, valores_eixo_x, valores_eixo_y, "k-");



