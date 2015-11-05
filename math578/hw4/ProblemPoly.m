%%
% Input data
n=10;
m=40;
randn('state',0);
u = linspace(-1,1,m);
v = 1./(5+40*u.^2) + 0.1*u.^3 + 0.01*randn(1,m);
A = vander(u');
A = A(:, m - n + [1:n]);

norms = [1, 1.5, 2, 2.5, 4, 8, inf];
l = length(norms);
max_residuals = zeros(1, l);

%%
% Construct polynomial for each norm and evaluate maximum norm of residual.

for i = 1 : length(norms)
    cvx_begin quiet
        variable x(n)
        minimize (norm(A*x - v', norms(i)))
    cvx_end
    vpol = x(1)*ones(1, length(u));
    for j = 2:n
      vpol = vpol.*u + x(j);
    end;
    max_residuals(i) = max(vpol - v);
end

%%
% Normalize and plot
max_residuals = max_residuals/max_residuals(l);

plot(norms(1:l-1), max_residuals(1:l-1), '-o')
xlabel('$p$')
ylabel('$\hat{r_p}$', 'interpreter', 'latex')
