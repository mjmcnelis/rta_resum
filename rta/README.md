**Gojko Vujanovic (c)**

This is my updated version of a [code](https://github.com/gvujan/Boltzmann_equation_solver_w_RTA_and_Bjorken_sym) provided by Gojko Vujanovic.

	L. Tinti, G. Vujanovic, J. Noronha, U. Heinz, Resummed hydrodynamic expansion for a plasma of particles interacting with fields, arXiv:1808.06436
	
The original [RTA code](http://personal.kent.edu/~mstrick6/code/) was originally developed by Michael Strickland

	W. Florkowski, R. Ryblewski, M. Strickland, Anisotropic Hydrodynamics for Rapidly Expanding Systems, arXiv:1304.0665.
	W. Florkowski, R. Ryblewski, M. Strickland, Testing viscous and anisotropic hydrodynamics in an exactly solvable case, arXiv:1305.7234.
	
You can use this code to evolve the exact conformal RTA Bjorken solution and also generate data for the figures in my thesis (Ch. 8)

	M. McNelis, arXiv:2105.06007

To compile and run, type make. The data is stored in `results`. To transfers the results to `bjorken_plot` for plotting, do

	sh copy_results_to_bjorken_plot.sh generator		# or ffe_1, ffe_2, attractor

You can edit the `PARAMETERS` macro in `rta.cpp` to generate data for each figure in the paper (and thesis).

