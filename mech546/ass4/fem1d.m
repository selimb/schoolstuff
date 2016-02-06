clear all;
%% Read Mesh and Problem Data
load('mesh.mat', 'nodeLocs', 'connectivity');
load('problem_data.mat', 'a', 'c', 'f');
numNodes = length(nodeLocs);
numElem = length(connectivity);

%% Construct global stiffness matrix and source terms
U = zeros(numNodes, 1);
K = zeros(numNodes, numNodes);
F = zeros(numNodes, 1);

for elem = 1:numElem
    conn = connectivity(elem,:);
    % Compute local stiffness and source
    leftNode = nodeLocs(conn(1),:);
    rightNode = nodeLocs(conn(2),:);
    h = rightNode(1) - leftNode(1);

    localK = (a/h)*[1 -1; -1 1] + (c*h/6)*[2 1; 1 2];
    localF = 0.5*f*h*[1; 1];  % Not handling Q for now

    % Add to global matrix
    K(conn, conn) = K(conn, conn) + localK;
    F(conn) = F(conn) + localF;
end
disp('Global Stiffness Matrix');
disp('');
disp(K);

%% Apply Boundary Conditions
% Hard-coded for now
Kc = K(2:numNodes-1, 2:numNodes-1);
Fc = F(2:numNodes-1);

%% Solve Linear System
% Hard-coded for now
U(2:numNodes-1) = Kc\Fc;
