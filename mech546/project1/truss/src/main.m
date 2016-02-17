clear all;
% Number of spatial dimensions is fixed.

%% Read Mesh and Problem Data
% Mesh
node_locs = csvread('node_locs.csv', 1, 0);
connectivity = csvread('connectivity.csv', 1, 0);
ndims = length(node_locs(1,:));  % Number of spatial dimensions

% Material properties (hard-coded)
% matprops = struct('A', 1, 'E', 1, 'RHO', 2770);
matprops = struct('A', 3225.8 * 10^-6, 'E', 69 * 10^9, 'RHO', 2770);
%E = 69 * 10^9;
%A = 3225.8 * 10^-6;
%RHO = 2770;
E = 1;
A = 1;
EA = E*A;

% Boundary Conditions
%   Applied Loads
L = csvread('loads.csv', 1, 0);
num_loads = size(L, 1);
loads = [];
for i = 1:num_loads
    loads = [loads, struct('node', L(i,1), 'dim', L(i,2), 'val', L(i,3))];
end

%   Essential BCs (homogeneous only)
fixed_nodes = csvread('fixed_nodes.csv', 1, 0)';

% Calculate numbers of stuff
num_nodes = length(node_locs);
num_elems = length(connectivity);
num_dofs = ndims*num_nodes;

%% Construct Stiffness Matrix
U = zeros(num_dofs, 1);
K = zeros(num_dofs, num_dofs);
F = zeros(num_dofs, 1);

for elem = 1:num_elems
    nodes = connectivity(elem,:);
    % Compute element stiffness
    [ Klocal, T ] = mk_stiff( node_locs(nodes,:), ndims );
    Kelem = EA*T'*Klocal*T;

    % Scatter into global
    sctr = mk_sctr(nodes, ndims);
    K(sctr, sctr) = K(sctr, sctr) + Kelem;
end
disp('Global Stiffness Matrix');
disp('');
K

%% Boundary Conditions
% Applied Loads
for load = loads
    sctr = mk_sctr(load.node, ndims);
    sctr = sctr(load.dim);
    F(sctr) = load.val;
end

% Essential Boundary Conditions
% Instead of deleting rows, we construct a mask. Beauty.
free_dofs = true(num_dofs,1);
for fixed_node = fixed_nodes
    sctr = mk_sctr(fixed_node, ndims);
    free_dofs(sctr) = false;
end
Kc = K(free_dofs, free_dofs);
Fc = F(free_dofs);

%% Solve Linear System
% Fixed DOFs are already 0.
disp('Condensed System')
disp('')
Kc
Fc
U(free_dofs) = Kc\Fc;
disp('Solve Complete')
U

%% Post-Computation
% Calculate Stresses
S = -1*ones(num_elems, 1);
for elem = 1:num_elems
    nodes = connectivity(elem,:);
    [ Kelem, T ] = mk_stiff(node_locs(nodes,:), ndims );
    sctr = mk_sctr(nodes, ndims);
    Ulocal = T*U(sctr);
    P = EA*Kelem*Ulocal;
    S(elem) = P(2)/A;
end
S
