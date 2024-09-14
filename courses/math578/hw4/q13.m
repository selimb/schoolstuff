function q13()
    X = [0 1; 0 -1; 0 0];
    y = [1 sqrt(3) 2 5];
    for i = 1 : length(y)
        X(3, 1) = y(i);
        go(X);
    end
end

function go(X)
    n = length(X);
    X = X';
    PTS = ones(2, 3);
    % 1-norm
    cvx_begin quiet
        variable a(2)
        minimize ( sum(dist(a, X, n)) )
    cvx_end
    PTS(:, 1) = a;
    % 2-norm
    cvx_begin quiet
        variable a(2)
        minimize ( sum( square_pos(dist(a, X, n)) ) )
    cvx_end
    PTS(:, 2) = a;
    % inf-norm
    cvx_begin quiet
        variable a(2)
        minimize ( max(dist(a, X, n)) )
    cvx_end
    PTS(:, 3) = a;
    styles = ['rs'; 'gd'; 'bx'];
    styles = cellstr(styles);  % Yay MATLAB...

    figure
    hold on
    for i = 1 : 3
        plot(PTS(1, i), PTS(2, i), styles{i});
    end
    legend('1-norm', '2-norm', 'inf-norm')
    xlim([0, 5])
    ylim([-1, 1])
    grid('off')
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    plot(X(1, :), X(2, :), 'ok')
end

function r = dist(a, X, n)
    r = norms(a*ones(1, n) - X, 2);
end
