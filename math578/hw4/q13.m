% TODO
% - Use CVX
% - Generate data
% - 2-3 functions with different data. Plot them all at the same time. Fuck it.
function q13main(X)
    X = [0 0; 0 1; 0 -1; 10 0];
    x = zeros(3, 2);

    x = [2 0; 2 1; 2 -1];

    norms = [1, 2, inf];
    styles = ['rs'; 'gd'; 'bx'];
    styles = cellstr(styles);  % Yay MATLAB...
    labels = [];

    figure
    hold on
    for i = 1 : 3
        plot(x(i, 1), x(i, 2), styles{i});
    end
    legend('1-norm', '2-norm', 'inf-norm')
    grid('off')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    plot(X(:, 1), X(:, 2), 'ok')
end
