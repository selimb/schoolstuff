% Generate Mesh
L = 1;
numElem = 5;
numNodes = numElem + 1;
x = linspace(0, L, numNodes);
y = zeros(1, length(x));
nodeLocs = [x; y]';
connectivity = [1:numNodes - 1; 2:numNodes]';
save('mesh.mat', 'nodeLocs', 'connectivity');

% Problem Data
a = -1;
c = 2;
f = 1;
save('problem_data.mat', 'a', 'c', 'f');

clear all;