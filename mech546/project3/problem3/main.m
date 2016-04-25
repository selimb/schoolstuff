clear all;
clc
%% Problem data
% Material properties
t = 0.1;
ANU = 0.1;
E = 2*10^6;

C = zeros(3);
C(1,1) = E/(1 - ANU^2);
C(1,2) = ANU*E/(1 - ANU^2);
C(3,3) = E/(2*(1 + ANU));
C(2,2) = C(1,1);
C(2,1) = C(1,2);

% Mesh
ndims = 2;
node_locs = csvread('node_locs.csv', 1, 0);
conn = csvread('connectivity.csv', 1, 0);

num_nodes = size(node_locs, 1);
num_elems = size(conn, 1);
num_dofs = ndims*num_nodes;

% Boundary conditions
%   Fixed boundaries
fixed_nodes = [ 3 4 5 5 6 7 ];
fixed_dims = [ 1 1 1 2 2 2 ];
num_fixed = length(fixed_nodes);
fixed_bcs = [];
for i = 1:num_fixed
    bc = struct('node', fixed_nodes(i), 'dim', fixed_dims(i));
    fixed_bcs = [fixed_bcs, bc];
end

%   Applied loads
p0 = 200;
l = 20;
f = 0.5*t*p0*l;
load_nodes = [ 1 2 3];
load_dims = [ 2 2 2];
load_vals = [ f 2*f f];
num_loads = length(load_nodes);
loads = [];
for i = 1:num_loads
    load = struct('node', load_nodes(i), 'dim', load_dims(i), 'val', load_vals(i));
    loads = [loads, load];
end

%% Construct Stiffness MatrixKc
K = zeros(num_dofs, num_dofs);
U = zeros(num_dofs, 1);
F = zeros(num_dofs, 1);
for elem = 1:num_elems
    nodes = conn(elem,:);
    % Compute element stiffness
    [A, B] = mk_B(node_locs(nodes,:));
    Kelem = t*A*(B'*C*B);

    % Scatter into global
    sctr = mk_sctr(nodes);
    K(sctr, sctr) = K(sctr, sctr) + Kelem;
end
disp('Global Stiffness Matrix');
K

%% Apply Boundary Conditions
% Applied Loads 
for load = loads
    sctr = mk_sctr(load.node, load.dim);
    F(sctr) = load.val;
end
disp('Global Source Term');
F

% Essential Boundary Conditions
% Instead of deleting rows, we construct a mask. Beauty.
free_dofs = true(num_dofs, 1);
for fixed = fixed_bcs
    sctr = mk_sctr(fixed.node, fixed.dim);
    free_dofs(sctr) = false;
end
Kc = K(free_dofs, free_dofs);
Fc = F(free_dofs);

%% Solve Linear System
disp('Condensed System')
Kc
Fc
% Fixed DOFs are already 0, only need to set the values of the free DOFs.
U(free_dofs) = Kc\Fc;
disp('Solve Complete')
U

% Write nodal displacements
Unodal = ones(num_nodes, ndims);
for node = 1:num_nodes
    sctr = mk_sctr(node);
    Unodal(node, :) = U(sctr);
end
csvwrite('Unodal.out', Unodal)

% Write elemental displacements
Uelem = ones(num_elems, ndims);
for elem = 1:num_elems
    nodes = conn(elem, :);
    nn = length(nodes);
    for dim = 1:ndims
        sctr = mk_sctr(nodes, dim);
        Uelem(elem, dim) = sum(U(sctr))/nn;
    end
end
csvwrite('Uelem.out', Uelem)

%% Calculate Stresses
S = zeros(num_elems, 3);
for elem = 1:num_elems
    nodes = conn(elem,:);
    [A, B] = mk_B(node_locs(nodes,:));
    sctr = mk_sctr(nodes);
    dof = U(sctr);
    S(elem, :) = C*B*dof;
end
disp('Stresses')
S
csvwrite('S.out', S)