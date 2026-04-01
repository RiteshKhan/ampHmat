# Adaptive mixed-precision hybrid hierarchical matrices

Hybrid hierarchical matrices are based on a *hybrid* admissibility condition: the standard admissibility condition is used at the coarser level of the hierarchy, and the weak admissibility condition is used at the finer level of the hierarchy.


## Please follow these simple steps to run the code.
Step 1.
```
git clone https://github.com/riteshkhan/ampHmat.git
```
Step 2.
```
cd ampHmat/Hmat
```

## The parameters in the main_file.m
1) d_dim = The dimension of the space. The tree will be generated accordingly. For example, if d_dim=3, the resulting tree will be an oct-tree.
1) nParticlesInLeafAlong1D = The maximum number of particles along 1D at the leaf level. The maximum number of particles at the leaf level N_max = nParticlesInLeafAlong1D^d_dim (this decides the depth of the tree $L$).  
3) LR_epsilon (ε) = Target accuracy in the low-rank compression.
4) kernel_choice = Choice of your desired kernel function (see kernel_menu.md file and select the number from there).
5) hodlr_level = This decides the switching level $\ell$ (where the admissibility conditions will be switched) in our hybrid admissibility condition. For example, if hodlr_level = 1, then switching level $\ell = L-1$.
6) eta = The parameter used in the standard admissibility condition. A different choice of eta will generate a different $\mathcal{H}$-matrix.
7) is_sym = The kernel matrix is symmetric or not (true/false).
8) u = Working precision for mixed precision H-matrix representation storage. 
9) u_mvp = Working precision for matrix-vector product with the resulting H-matrix representation.
10) The code can be run in mixed precision as well as uniform precision settings (more info can be found in the main_file.m).



<p align="center">
  <img src="Hmat/Images/hmat1.jpg" alt="H-matrix example" width="1000"><br>
  <em>
    <em>Figure 1: Colors indicate the precision used in each block for target accuracy ε = 10⁻⁴</em>.
  </em>
</p>

<br>

<p align="center">
  <img src="Hmat/Images/hmat2.jpg" alt="H-matrix example" width="1000"><br>
  <em>
    <em>Figure 2: Colors indicate the precision used in each block for target accuracy ε = 10⁻²</em>.
  </em>
</p>

Thus, we can see that almost all the blocks can be stored in a precision lower than the working precision. Notably, this does not come at the expense of accuracy.

## Acknowledgement 
This work is supported by the European Union (ERC, inEXASCALE, 101075632) and the Charles University Research Centre program No. UNCE/24/SCI/005. 
