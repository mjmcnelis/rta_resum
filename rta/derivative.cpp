
void compute_first_derivative(double *fprime, double *f, int Nx, double dx)
{
	for(int i = 0; i < Nx; i++)
  {
    if(i == 0)
    {
      fprime[i] = (f[i+1] - f[i]) / dx;          // forward difference
    }
    else if(i == Nx - 1)
    {
      fprime[i] = (f[i] - f[i-1]) / dx;          // backward difference
    }
    else
    {
      fprime[i] = (f[i+1] - f[i-1]) / (2.*dx);   // central difference
    }
  }
}
