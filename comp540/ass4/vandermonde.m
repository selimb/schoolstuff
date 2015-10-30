function V = vander(m, n)
    V = zeros(m, n);
    for i = 1 : m
        for j = 1 : n
            V(i, j) = (j/n)^(i - 1);
        end
    end
end
