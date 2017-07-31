# distance_sparcifier

### Description
distance_sparcifier is a C++ code to sparcify distance matrices as a preprocessing step for the computation of sparse Vietoris-Rips persistence barcodes. 
The code for reading distance matrices is based on [Ripser](https://github.com/Ripser). 

### Building

distance_sparcifier requires a C++11 compiler. Here is how to obtain and build distance_sparcifier:

```sh
git clone https://github.com/blasern/distance_sparcifier
cd distance_sparcifier
make
```

### Use

Here is how to run distance_sparcifier:

```sh
./distance_sparcifier examples/sphere_3_192.lower_distance_matrix 
```
