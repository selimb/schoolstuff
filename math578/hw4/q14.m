%% Read data and classify points
load('discrim.mat');
U = Y(Y(:, 3) == 0, 1:2)';
V = Y(Y(:, 3) == 1, 1:2)';

%% CVX
cvx_begin quiet
    variables c(2) b(1) t(1)
    maximize (t)
    c'*U - b >= t;
    c'*V - b <= -t;
%   norm(c) <= 1;
cvx_end

%% Plot
xc = linspace(0, 1);
% x*c1 + y*c2 = b
yc = (-c(1)*xc + b)/c(2);
figure
hold on
plot(U(1, :), U(2, :), 'ro')
plot(V(1, :), V(2, :), 'bo')
plot(xc, yc, 'k--')

%% Print
t
c
b
