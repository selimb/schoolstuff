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
% I'm personally not a fan of this mapping because it makes the scattering
% less obvious, but whatever.
%
% Consequently, the `i` component of displacement at node `I` is located at:
%   ndims*(I - 1) + i
%
% Example:
% >> mk_sctr([1 3], 2)
%   [1 2 5 6]
%
% >> mk_sctr([1 3], 3)
%   [1 2 3 7 8 9]
function [ sctr ] = mk_sctr(nodeIDs, ndims)
    % Pre-allocating would be faster, but requires non-trivial
    % indexing calculations or, even worse, an additional loop.
    sctr = [];
    for I = nodeIDs
        dofs = ndims*(I - 1) + [1:ndims];
        sctr = [sctr dofs];
    end
end
