#!/usr/bin/env sh
echo 'Enter directory name:'
read dir

mkdir -p $dir
echo Connectivities > $dir/connectivity.csv
echo "node" > $dir/fixed_nodes.csv
echo "node, dim, val" > $dir/loads.csv
echo "x, y, z" > $dir/node_locs.csv
