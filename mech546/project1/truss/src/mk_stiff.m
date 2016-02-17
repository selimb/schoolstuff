function [ K ] = mk_stiff(node_locs, ndims)
% Calculate and construct "specific" element stiffness matrix.
%   "specific" refers to the fact that material properties are not
%   accounted for.
% `node_locs` is a (2 x ndims) array
    diff = node_locs(2,:) - node_locs(1,:);
    if ndims == 2
        K = stiffness_2d(diff);
    elseif ndims == 3
        K = stiffness_3d(diff);
    else
        error('Cannot calculate stiffness.')
    end
end

function [ K ] = stiffness_2d(diff)
    dx = diff(1);
    dy = diff(2);
    L = sqrt( sum(diff.^2) );
    c = dx/L;
    s = dy/L;
    c2 = c*c; cs = c*s; s2 = s*s;

    A = [c2 cs;
         cs s2];
    K = (1/L) * [ A -A;
                 -A  A];
end

function [ K ] = stiffness_3d(diff)
    dx = diff(1);
    dy = diff(2);
    dz = diff(3);
    L = sqrt( sum(diff.^2) );
    cx = dx/L;
    cy = dy/L;
    cz = dz/L;
    cx2 = cx*cx; cy2 = cy*cy; cz2=cz*cz;
    cxy = cx*cy; cxz = cx*cz; cyz = cy*cz;

    A = [cx2 cxy cxz;
         cxy cy2 cyz;
         cxz cyz cz2];
    K = (1/L) * [ A -A;
                 -A  A];
end
