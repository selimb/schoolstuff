function [A, B] = mk_B(node_locs)
    % `node_locs` is a (3x2) array
    alphas = zeros(3,1);
    betas = zeros(3,1);
    gammas = zeros(3,1);
    x = zeros(5,1); y = zeros(5,1);
    x(1:3) = node_locs(:,1); y(1:3) = node_locs(:,2);
    x(4) = x(1); x(5) = x(2); y(4) = y(1); y(5) = y(2);
    for i = 1:3
        j = i + 1;
        k = i + 2;
        alphas(i) = x(j)*y(k) - x(k)*y(j);
        betas(i) = y(j) - y(k);
        gammas(i) = -(x(j) - x(k));
    end
    A = 0.5*sum(alphas);
    B = zeros(3,6);
    B(1, 1:2:end) = betas;
    B(2, 2:2:end) = gammas;
    B(3, 1:2:end) = gammas;
    B(3, 2:2:end) = betas;
    B = B/(2*A);
end