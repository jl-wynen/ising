# Ising
This repository is a collection of a number of different implementations of the [Ising Model](https://en.wikipedia.org/wiki/Ising_model).
There are simple implementations of the basic 2D model without external magnetic field in various languages as well as one implementation of the N dimensional model in C++.
The former are meant to serve as a comparison of different programming languages with respect to ease of use and run efficiency.
See the sub directories [comparison](/comparison) and [n-dimensional](n-dimensional) for more information.

## Model
The Ising Model is given by the Hamiltonian (assuming no external magnetic field)  
<a href="https://www.codecogs.com/eqnedit.php?latex=H(s)&space;=&space;-J&space;\sum_{\langle&space;i,j&space;\rangle}\,&space;s_i&space;s_j&space;-&space;h&space;\sum_{i}\,&space;s_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H(s)&space;=&space;-J&space;\sum_{\langle&space;i,j&space;\rangle}\,&space;s_i&space;s_j&space;-&space;h&space;\sum_{i}\,&space;s_i" title="H(s) = -J \sum_{\langle i,j \rangle}\, s_i s_j - h \sum_{i}\, s_i" /></a>
where s is a configuration of spins with s_i=+1,-1. The angle brackets denote nearest neighbours.
The units are normalised such that k_B T = 1 which implies that J and h  are dimensionless.
