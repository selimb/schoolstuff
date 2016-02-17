% Test stiffness matrix constructions
iseq = @(A,B, tol) ( all(all(abs(A-B) <= tol)) );

%% 2-D
ndims = 2;

% From Example 6.4.1
tol = 1e-7;
K = mk_stiff([0 0; 1 0], ndims);
Kex = [1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];
assert (iseq(K, Kex, tol), 'elem1');

K = mk_stiff([0 0; 0 1], ndims);
Kex = [0 0 0 0; 0 1 0 -1; 0 0 0 0; 0 -1 0 1];
assert (iseq(K, Kex, tol), 'elem2');

K = mk_stiff([0 0; 1 1], ndims);
c = cosd(45)^2;
Kex = [c c -c -c; c c -c -c; -c -c c c; -c -c c c]/sqrt(2);
assert (iseq(K, Kex, tol), 'elem3');

% From http://www.ce.memphis.edu/7117/notes/presentations/chapter_03a.pdf
% p.39
K = mk_stiff([0 0; 3 4], ndims);
Aex = [0.36 0.48; 0.48 0.64];
Kex = [Aex -Aex; -Aex Aex]/5;
assert (iseq(K, Kex, 1e-4), 'elemx');
%% 3-D
ndims = 3;

% From http://www.ce.memphis.edu/7117/notes/presentations/chapter_03a.pdf
% p.46
K = mk_stiff([72 0 0; 0 36 0], 3);
Aex = [0.79 -0.4 0; -0.4 0.2 0; 0 0 0];
Kex = (1/80.5)*[Aex -Aex; -Aex Aex];
assert (iseq(K, Kex, 1e-3), '3d elem1');

K = mk_stiff([72 0 0; 0 36 72], 3);
Aex = [0.45 -0.22 -0.45; -0.22 0.11 0.45; -0.45 0.45 0.45];
Kex = (1/108.0)*[Aex -Aex; -Aex Aex];
assert (iseq(K, Kex, 1e-2), '3d elem2');