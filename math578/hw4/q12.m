%% Input data
M = ones(4, 1);
b = [5.; -5.; 1.; -1.];

norms = [2, inf];
num_norms = length(norms);
best_x = zeros(1, num_norms);

%% Construct polynomial for each norm and evaluate maximum norm of residual.
for i = 1 : num_norms
    cvx_begin quiet
        variable x
        minimize (norm(M*x - b, norms(i)))
    cvx_end
    best_x(i) = x;
end

%% Brute force
x = linspace(-5, 5, 101);
norm_results = zeros(num_norms, length(x));

for i = 1 : length(x)
    r = M*x(i) - b;
    for n = 1 : num_norms
        norm_results(n, i) = norm(r, norms(n));
    end
end

%% Output and plot
fprintf(1, 'The best x according to CVX are:\n')
disp(best_x)

legend_names = cell(1, num_norms);
figure
hold on
for n = 1 : num_norms
    plot(x, norm_results(n, :))
    legend_names{n} = sprintf('%s norm', num2str(norms(n)));
end
legend(legend_names, 'Location', 'North')
xlabel('$x$', 'interpreter', 'latex')
ylabel('$||r||$', 'interpreter', 'latex')
title('Variation of residual with $x$ for different norms', 'interpreter', 'latex')
