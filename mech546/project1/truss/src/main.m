clear all;
% Number of spatial dimensions is fixed.

%% Read Mesh and Problem Data
NDIMS = 3;  % Number of spatial dimensions
nodeLocs = csvread('nodeLocs.csv', 1, 0);
connectivity = csvread('connectivity.csv', 1, 0);

A = 3225.8 * 10^-6;
E = 69 * 10^9;
RHO = 2770;

% Boundary Conditions (hard-coded)
%   Load at node 1 in z-direction of 10 000N
loads = [ struct('node', 1, 'dim', 3, 'val', 10000) ];
%   Essential BCs
fixed_nodes = [7 8 9 10];

num_nodes = length(nodeLocs);
num_elems = length(connectivity);
num_dofs = NDIMS*num_nodes;

%% Construct Stiffness Matrix
U = zeros(num_dofs, 1);
K = zeros(num_dofs, num_dofs);
F = zeros(num_dofs, 1);

for elem = 1:num_elems
    nodes = connectivity(elem,:);
    % Compute element stiffness


    % Scatter into global
    sctr = mk_sctr(nodes, NDIMS);
    K(sctr, sctr) = K(sctr, sctr) + localK;
end
disp('Global Stiffness Matrix');
disp('');
disp(K);

%% Boundary Conditions
% Applied Loads
for load = loads
    sctr = mk_sctr(load.node, NDIMS);
    sctr = sctr(load.dim);
    F(sctr) = load.val;
end

% Essential Boundary Conditions
% Instead of deleting rows, we construct a mask. Beauty.
free_dofs = true(num_dofs,1);
for fixed_node = fixed_nodes
    sctr = mk_sctr(fixed_node, NDIMS);
    free_dofs(sctr) = false;
end
Kc = K(free_dofs, free_dofs);
Fc = F(free_dofs);

%% Solve Linear System
% Fixed DOFs are already 0.
U(free_dofs) = Kc\Fc;
