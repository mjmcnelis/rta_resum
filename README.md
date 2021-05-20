**RTA Resum (c) Mike McNelis** 

Created on 12/11/2019 by Mike McNelis \
Last edited on 5/19/2021 by Mike McNelis

# Summary

This repository contains supplemental material for my hydrodynamic generator paper and PhD thesis:

	M. McNelis, U. Heinz, Phys. Rev. C 101 054901 (2020)
	M. McNelis, arXiv:2105.06007

The directory `rta` contains a C++ code that evolves the exact solution of the conformal RTA Boltzmann equation (and several hydrodynamic approximations). I used this code to generate the plots in my thesis. 

The notebook in `generator_expansion` demonstrates the correspondence betweeen the hydrodynamic generator and RTA Chapman-Enskog expansion. A more concrete proof can be found in Appendix E.1 of the thesis.  

The notebook in `standard_df_corrections` computes the Chapman-Enskog df corrections up to third order. These are also used to evaluate the terms in the hydrodynamic generator expansion.

