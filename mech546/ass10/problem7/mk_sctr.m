% Construct scatter vector used for operator scattering (Assembly).
%
% Scattering is the step where we "plug in" a local operator, e.g.
% stiffness matrix, into the global one. The returned scatter vector
% `sctr` can be used as follows to scatter a local 1-D array:
%     Fglobal(sctr) = Fglobal(sctr) + Flocal
% and into a matrix:
%     Kglobal(sctr,sctr) = Kglobal(sctr,sctr) + Klocal
%
% This is entirely dependent on the chosen DOF mapping. To be consitent
% with the book and lecture notes, we first alternate between each
% spatial component and then alternate between nodes in ascending order:
%   U = [u1x u1y u1z u2x u2y u2z ... uN1 uNy uNz ]
%
% Consequently, the `i` component of displacement at node `I` is located at:
%   ndims*(I - 1) + i
%
% For example:
% >> mk_sctr([1 3])
%   [1 2 5 6]
%
% Moreover, if one wishes to obtain the scatter vector in a particular
% component, it can be passed into the `dim` argument:
%
% >> mk_sctr([1 3], 2)
%   [2 6]

function [ sctr ] = mk_sctr(nodeIDs, dim)
    do_dim = true;
    if nargin == 1
        do_dim = false;
    end
    % Pre-allocating would be faster, but requires non-trivial
    % indexing calculations or, even worse, an additional loop.
    ndims = 2;
    sctr = [];
    for I = nodeIDs
        dofs = ndims*(I - 1) + [1:ndims];
        if do_dim
            dofs = dofs(dim);
        end
        sctr = [sctr dofs];
    end
end