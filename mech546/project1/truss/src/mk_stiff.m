function [ Klocal, T ] = mk_stiff(node_locs, ndims)
% Calculate local K and transformation matrix T
%   "specific" refers to the fact that material properties are not
%   accounted for.
% `node_locs` is a (2 x ndims) array
    Klocal = [1 -1; -1 1];
    diff = node_locs(2,:) - node_locs(1,:);
    if ndims == 2
        [ L, T ] = stiffness_2d(diff);
    elseif ndims == 3
        [ L, T ] = stiffness_3d(diff);
    else
        error('Cannot calculate stiffness.')
    end
    Klocal = (1/L)*Klocal;
end

function [ L, T ] = stiffness_2d(diff)
    dx = diff(1);
    dy = diff(2);
    L = sqrt( sum(diff.^2) );
    c = dx/L;
    s = dy/L;
    Tq = [c s];
    z = zeros(size(Tq));
    T = [Tq z; z Tq];
    % Tq = [c s; -s c];
    % z = zeros(size(Tq));
    % T = [Tq z; z Tq];
end

function [ L, T ] = stiffness_3d(diff)
    dx = diff(1);
    dy = diff(2);
    dz = diff(3);
    L = sqrt( sum(diff.^2) );
    cx = dx/L;
    cy = dy/L;
    cz = dz/L;
    Tq = [cx cy cz];
    z = zeros(size(Tq));
    T = [Tq z; z Tq]
end

