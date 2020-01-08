
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


void compute_second_derivative(double *f2prime, double *f, int Nx, double dx)
{
  for(int i = 0; i < Nx; i++)
  {
    if(i == 0)
    {
      f2prime[i] = (f[i+2] - 2.*f[i+1] + f[i]) / (dx*dx);    // forward difference
    }
    else if(i == Nx - 1)
    {
      f2prime[i] = (f[i] - 2.*f[i-1] + f[i-2]) / (dx*dx);    // backward difference
    }
    else
    {
      f2prime[i] = (f[i+1] - 2.*f[i] + f[i-1]) / (dx*dx);    // central difference
    }
  }
}

void compute_third_derivative(double *f3prime, double *f, int Nx, double dx)
{
  for(int i = 0; i < Nx; i++)
  {
    if(i <= 2)
    {
      f3prime[i] = (f[i+3] - 3.*f[i+2] + 3.*f[i+1] - f[i]) / (dx*dx*dx);      // forward difference
    }
    else if(i >= Nx - 3)
    {
      f3prime[i] = (f[i] - 3.*f[i-1] + 3.*f[i-2] - f[i-3]) / (dx*dx*dx);      // backward difference
    }
    else
    {
      f3prime[i] = (f[i+3] - 3.*f[i+1] + 3.*f[i-1] - f[i-3]) / (8.*dx*dx*dx); // central difference
    }
  }
}


void compute_fourth_derivative(double *f4prime, double *f, int Nx, double dx)
{
  for(int i = 0; i < Nx; i++)
  {
    if(i <= 1)
    {
      f4prime[i] = (f[i+4] - 4.*f[i+3] + 6.*f[i+2] - 4.*f[i+1] + f[i]) / (dx*dx*dx*dx);   // forward difference
    }
    else if(i >= Nx - 2)
    {
      f4prime[i] = (f[i] - 4.*f[i-1] + 6.*f[i-2] - 4.*f[i-3] + f[i-4]) / (dx*dx*dx*dx);   // backward difference
    }
    else
    {
      f4prime[i] = (f[i+2] - 4.*f[i+1] + 6.*f[i] - 4.*f[i-1] + f[i-2]) / (dx*dx*dx*dx);   // central difference
    }
  }
}







