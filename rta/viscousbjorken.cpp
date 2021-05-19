#include <math.h>
#include <cmath>
#include <iostream>


double de_dt(double e, double pi, double t, double prefactor, double etas, double third)
{
  return (pi  -  4./3. * e) / t;
}


double dpi_dt(double e, double pi, double t, double prefactor, double etas, double third)
{
  double p = e/3.;
  double T = pow(e / prefactor, 0.25);
  double taupi_inv = T / (5. * etas);

  return - taupi_inv * pi  +  16. * p / (15. * t)  -  38. * pi / (21. * t)  -  third * 18. * pi * pi / (49. * p * t);
}


void run_second_order_viscous_bjorken(double T0, double pl0, double t0, double dt, int Nt, double prefactor, double etas)
{
  double t  = t0;
  double e = prefactor * pow(T0, 4);
  double pi = e/3.  -  pl0;

  double third = 0;

  FILE *shear;
  shear = fopen("results/pibar_DNMR.dat", "w");

  for(int i = 1; i < Nt; i++)
  {
    fprintf(shear, "%.6e\t%.6e\n", t, 3. * pi / e);     // pi / peq

    double e1  = dt *  de_dt(e, pi, t, prefactor, etas, third);
    double pi1 = dt * dpi_dt(e, pi, t, prefactor, etas, third);

    double e2  = dt *  de_dt(e + e1/2., pi + pi1/2., t + dt/2., prefactor, etas, third);
    double pi2 = dt * dpi_dt(e + e1/2., pi + pi1/2., t + dt/2., prefactor, etas, third);

    double e3  = dt *  de_dt(e + e2/2., pi + pi2/2., t + dt/2., prefactor, etas, third);
    double pi3 = dt * dpi_dt(e + e2/2., pi + pi2/2., t + dt/2., prefactor, etas, third);

    double e4  = dt *  de_dt(e + e3, pi + pi3, t + dt, prefactor, etas, third);
    double pi4 = dt * dpi_dt(e + e3, pi + pi3, t + dt, prefactor, etas, third);

    e  += (e1   +  2. * e2   +  2. * e3   +  e4) / 6.;
    pi += (pi1  +  2. * pi2  +  2. * pi3  +  pi4) / 6.;

    t += dt;
  }

  fclose(shear);
}


void run_third_order_viscous_bjorken(double T0, double pl0, double t0, double dt, int Nt, double prefactor, double etas)
{
  double t  = t0;
  double e = prefactor * pow(T0, 4);
  double pi = e/3.  -  pl0;

  double third = 1;

  FILE *shear;
  shear = fopen("results/pibar_Jaiswal.dat", "w");

  for(int i = 1; i < Nt; i++)
  {
    fprintf(shear, "%.6e\t%.6e\n", t, 3. * pi / e);     // pi / peq

    double e1  = dt *  de_dt(e, pi, t, prefactor, etas, third);
    double pi1 = dt * dpi_dt(e, pi, t, prefactor, etas, third);

    double e2  = dt *  de_dt(e + e1/2., pi + pi1/2., t + dt/2., prefactor, etas, third);
    double pi2 = dt * dpi_dt(e + e1/2., pi + pi1/2., t + dt/2., prefactor, etas, third);

    double e3  = dt *  de_dt(e + e2/2., pi + pi2/2., t + dt/2., prefactor, etas, third);
    double pi3 = dt * dpi_dt(e + e2/2., pi + pi2/2., t + dt/2., prefactor, etas, third);

    double e4  = dt *  de_dt(e + e3, pi + pi3, t + dt, prefactor, etas, third);
    double pi4 = dt * dpi_dt(e + e3, pi + pi3, t + dt, prefactor, etas, third);

    e  += (e1   +  2. * e2   +  2. * e3   +  e4) / 6.;
    pi += (pi1  +  2. * pi2  +  2. * pi3  +  pi4) / 6.;

    t += dt;
  }

  fclose(shear);
}


