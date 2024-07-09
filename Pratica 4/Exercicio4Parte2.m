clear
A = rand(3)

C = A'*A;

[V, sigma_quad] = eig(C);

sigma = sqrt(diag(sigma_quad));

Sigma = zeros(size(A));
for i = 1:length(sigma)
    Sigma(i, i) =sigma(i);
end

U = A * V * inv(Sigma)
Sigma
V'


U*Sigma*V'

