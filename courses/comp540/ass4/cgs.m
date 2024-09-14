function [Q, R] = cgs(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);

    for k = 1 : n
        for i = 1 : k - 1
            R(i, k) = Q(:, i)'*A(:, k);
        end
        for i = 1 : k - 1
            A(:, k) = A(:, k) - R(i, k).*Q(:, i);
        end
        R(k, k) = norm(A(:, k));
        Q(:, k) = A(:, k)/R(k, k);
    end
end
