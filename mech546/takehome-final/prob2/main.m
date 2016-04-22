clc
%% Problem data
% Material properties
h = 1;
v = 0.3;
E = 30*10^6;

C = zeros(3);
C(1,1) = E/(1 - v^2);
C(1,2) = v*E/(1 - v^2);
C(3,3) = E/(2*(1 + v));
C(2,2) = C(1,1);
C(2,1) = C(1,2);

% Mesh
ndims = 2;
node_locs = [ 0 0; 0 10; 20 10; 20 0];
connectivity = [1 4 3; 3 2 1];

num_nodes = size(node_locs, 1);
num_elems = size(connectivity, 1);
num_dofs = ndims*num_nodes;

% Boundary conditions
%   Fixed boundaries
fixed_nodes = [ 1 1 2 2];
fixed_dirs = [ 1 2 1 2];
fixed_bcs = [];
for i = 1:length(fixed_nodes)
    bc = struct('node', fixed_nodes(i), 'dim', fixed_dirs(i));
    fixed_bcs = [fixed_bcs, bc];
end


%% Construct Stiffness Matrix
K = zeros(num_dofs, num_dofs);
U = zeros(num_dofs, 1);
for elem = 1:num_elems
    nodes = connectivity(elem,:);
    % Compute element stiffness
    [A, B] = mk_B(node_locs(nodes,:));
    Kelem = h*A*(B'*C*B);

    % Scatter into global
    sctr = mk_sctr(nodes);
    K(sctr, sctr) = K(sctr, sctr) + Kelem;
end
K;

%% Apply Boundary Conditions
% Applied Loads (hard-coded)
p0 = 1000;
b = 10;
h = 1;
f = p0*b*h/2;
Fc = [f; 0; f; 0];

% Essential Boundary Conditions
% Instead of deleting rows, we construct a mask. Beauty.
free_dofs = true(num_dofs,1);
for fixed = fixed_bcs
    sctr = mk_sctr(fixed.node);
    sctr = sctr(fixed.dim);
    free_dofs(sctr) = false;
end
Kc = K(free_dofs, free_dofs);
Kc

%% Solve Linear System
disp('Condensed System')
Kc
Fc
% Fixed DOFs are already 0, only need to set the values of the free DOFs.
U(free_dofs) = Kc\Fc;
disp('Solve Complete')
U

%% Calculate Stresses
for elem = 1:num_elems
    nodes = connectivity(elem,:);
    [A, B] = mk_B(node_locs(nodes,:));
    sctr = mk_sctr(nodes);
    dof = U(sctr);
    s = C*B*dof
end