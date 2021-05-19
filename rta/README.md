Gojko Vujanovic (c) 2018

The code solves the conformal RTA Bjorken solution. It was originally provided by Gojko Vujanovic (see https://github.com/gvujan/Boltzmann_equation_solver_w_RTA_and_Bjorken_sym)

	L. Tinti, G. Vujanovic, J. Noronha, U. Heinz, Resummed hydrodynamic expansion for a plasma of particles interacting with fields, arXiv:1808.06436

The code in the link is an improved version of a code originally developed by Michael Strickland

	W. Florkowski, R. Ryblewski, M. Strickland, Anisotropic Hydrodynamics for Rapidly Expanding Systems, arXiv:1304.0665.
	W. Florkowski, R. Ryblewski, M. Strickland, Testing viscous and anisotropic hydrodynamics in an exactly solvable case, arXiv:1305.7234.

To compile and run, type make. The results are stored in `results`. To transfers the results to `bjorken_plot` for plotting, do

	sh copy_results_to_bjorken_plot generator		# or ffe_1, ffe_2, attractor

You can edit the `PARAMETERS` macro in `rta.cpp` to generate data for each figure in the paper (and thesis).

