# distance_sparcifier

### Description
distance_sparcifier is a C++ code to sparcify distance matrices as a preprocessing step for the computation of sparse Vietoris-Rips persistence barcodes. The code for reading distance matrices, the help and the examples are based on [Ripser](https://github.com/Ripser). 

Sparcification is based on [Cavanna et al](https://arxiv.org/abs/1506.03797). Let P be a furthest point sampling of a point cloud with insertion radii λi. 
Given δ with 0 < δ < 1 we define r0 = r1 = ∞ and ri = 2λ_{i−1} /(1 − δ ) for i = 2, ... , n.
In order to sparsify the edge list for the full Rips complex we set the distance between every edge (i, j) with d(i, j) ≥ min(ri , rj) to the maximal distance in P.


### Building

distance_sparcifier requires a C++11 compiler. Here is how to obtain and build distance_sparcifier:

```sh
git clone https://github.com/blasern/distance_sparcifier
cd distance_sparcifier
make
```

### Use

Here is how to run distance_sparcifier and writing the output to a lower distance matrix file:

```sh
./distance_sparcifier --interleaving 0.5 examples/sphere_3_192.lower_distance_matrix > examples/sphere_3_192_sparse_0.5.lower_distance_matrix
./distance_sparcifier --interleaving 0.3 examples/sphere_3_192.lower_distance_matrix > examples/sphere_3_192_sparse_0.3.lower_distance_matrix
./distance_sparcifier --interleaving 0.1 examples/sphere_3_192.lower_distance_matrix > examples/sphere_3_192_sparse_0.1.lower_distance_matrix
```

After this preprocessing step, run ripser as usual: 
```sh
./ripser examples/sphere_3_192.lower_distance_matrix
./ripser examples/sphere_3_192_sparse_0.5.lower_distance_matrix 
./ripser examples/sphere_3_192_sparse_0.3.lower_distance_matrix
./ripser examples/sphere_3_192_sparse_0.1.lower_distance_matrix
```
