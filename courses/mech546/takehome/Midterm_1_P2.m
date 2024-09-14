%% Setup
% Constants
a0 = 180*10^6;
f0 = 1000;
L = 10;

% Initialize
C = zeros(6,6);

%% c)
for N = 1:6
    % Construct A and b
    A = zeros(N,N);
    b = zeros(N,1);
    for i = 1:N
        b(i) = f0*L^(i+1)/(i+1);
        for j = 1:N
            A(i,j) = a0*i*j*L^(i+j-1) * (2/(i+j-1) - 1/(i+j));
        end
    end
    % Solve for coefficients
    coeffs = A\b;
    C(N,1:N) = coeffs;
end

C

%% d)
hold on;
x = linspace(0, L);
for N = 1:6
    % Construct polynomial
    U = zeros(1,length(x));
    coeffs = C(N,1:N);
    for i = 1:N
        U = U + coeffs(i)*(x.^i);
    end
    plot(x, U);
end
xlabel('x')
ylabel('u')
legend('1', '2', '3', '4', '5', '6')