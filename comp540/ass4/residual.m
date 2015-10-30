function s = residual(A, Q, R)
    s = sprintf('Residual: %.1e', norm(A - Q*R)/norm(A));
end
