function s = orthoerr(Q1)
    [m, n] = size(Q1);
    I = eye(n);
    s = sprintf('Orthoerr: %.1e', norm(Q1'*Q1 - I));
end
