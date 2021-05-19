#include <math.h>
#include <cmath>
#include <iostream>
#include "include/hypergeometric.h"


double compute_A(double x, double L04)
{
  double dx = 0.01;

  double HL = 0;
  double H240 = 0;

  if(fabs(x - 1.) <= dx)                    // taylor expansion around x = 1
  {
    HL = 0.004372373412624753 + x*(-0.04433956136745337 + x*(0.2078416939099339 +
         x*(0.9680100429323358 + x*(-0.7695771720535178 + x*
                (0.42777119681106374 + x*(-0.15963396768329577 + (0.03589393063070768 - 0.0036718699257309944*x)*x))))));


    H240 = -0.005087625772544818 + x*(0.05695884534842821 + x*(-0.3018911943816236 +
         x*(1.019465776582752 + x*(-0.5047550486743451 + x*(0.16392338295963027 +
                  x*(-0.030830985240682104 + (0.0020622455600826832 + 0.0001546036183062223*x)*x))))));
  }
  else if(x > 0. && x < 1. - dx)            // 0 < x < 1
  {
    double x2 = x * x;
    double x4 = x2 * x2;
    double z = sqrt(1./x2  -  1.);
    double t = atan(z) / z;

    HL = (-x2 + t) * x2 / (1. - x2);
    H240 = (2.  +  x2  -  3. * t) * x4 / (x2 - 1.) / (x2 - 1.);
  }
  else if(x > 1. + dx)                     // x > 1
  {
    double x2 = x * x;
    double x4 = x2 * x2;
    double z = sqrt(1. - 1./x2);
    double t = atanh(z) / z;

    HL = (-x2 + t) * x2 / (1. - x2);
    H240 = (2.  +  x2  -  3. * t) * x4 / (x2 - 1.) / (x2 - 1.);
  }
  else
  {
    printf("compute_A error: x = %lf is out of bounds", x);
    exit(-1);
  }

  return L04 / (M_PI * M_PI) * (4. * HL  -  1.5 * H240);
}


double compute_I240(double pl, double e)
{
  double daL = 0.01;
  double x   = pl / e;

  double x2  = x   * x;
  double x3  = x2  * x;
  double x4  = x3  * x;
  double x5  = x4  * x;
  double x6  = x5  * x;
  double x7  = x6  * x;
  double x8  = x7  * x;
  double x9  = x8  * x;
  double x10 = x9  * x;
  double x11 = x10 * x;
  double x12 = x11 * x;
  double x13 = x12 * x;
  double x14 = x13 * x;
  double x15 = x14 * x;
  double x16 = x15 * x;
  double x17 = x16 * x;
  double x18 = x17 * x;
  double x19 = x18 * x;
  double x20 = x19 * x;             // x = pl / e
  double x21 = x20 * x;             // z = xi (from aniso hydro)
  double x22 = x21 * x;             // aL = 1/sqrt(1+xi)

  double aL, aL2, xi;

  if(x > 1./3.)
  {
    xi = (33920.75424130814 - 72210.38155086128*x - 58150.373221605056*x2 - 123258.68865968155*x3 + 14269.18945991164*x4 + 170364.85584208343*x5 +
      169910.50665817957*x6 + 64799.56668758523*x7 + 262850.50962558796*x8 + 35449.81323782106*x9 - 248808.80620651352*x10 -
      288836.0950617432*x11 - 55525.59817904083*x12 - 249947.48251438234*x13 - 310420.8253593438*x14 - 48317.55700989516*x15 +
      522691.7302032236*x16 + 527504.8150488662*x17 + 219759.88337782127*x18 - 187603.57642353655*x19 - 199506.45061878706*x20 -
      611077.848257917*x21 + 432142.0012199023*x22)/
      (-192.35667843131404 + 43491.082537769165*x + 83354.47448892899*x2 + 45103.07343085356*x3 - 105414.36804542418*x4 + 140186.71296754244*x5 -
      531082.9994828509*x6 - 85658.91194364589*x7 + 377783.60198413196*x8 - 339045.0410056553*x9 - 95837.02795785779*x10 +
      284537.5663725089*x11 + 703062.1998023012*x12 + 223019.9316692852*x13 - 501784.5491947427*x14 - 145230.1534789184*x15 +
      55948.62853147295*x16 + 49679.34805386173*x17 - 641771.3022609851*x18 - 274804.41532698454*x19 + 726388.8998660464*x20 +
      350014.57800287893*x21 - 361748.9148710701*x22);

    if(!(xi >= -0.99999999 && xi <= 1.e-8))
    {
      xi = fmax(-0.99999999, fmin(xi, 1.e-8));
    }
    aL = 1. / sqrt(1. + xi);
    aL2 = 1. / (1. + xi);
  }
  else
  {
    aL = (2.372796737893896e-62 + 9.355496760751141e-54*x + 3.15985529218801e-46*x2 + 2.29804656071578e-39*x3 + 4.8654069671748624e-33*x4 +
      3.4686835009134695e-27*x5 + 9.052410236842743e-22*x6 + 9.132309729581051e-17*x7 + 3.705485165853083e-12*x8 + 6.240802836268058e-8*x9 +
      0.00044799689605487286*x10 + 1.4025011569370325*x11 + 1953.0793537979494*x12 + 1.229812256787706e6*x13 + 3.543561225712354e8*x14 +
      4.697330865356272e10*x15 + 2.8499566740003765e12*x16 + 7.731610782177606e13*x17 + 8.841791912264315e14*x18 + 3.673281425421166e15*x19 +
      3.042059896930142e15*x20 - 3.4368817938638095e15*x21 - 8.169907788507815e14*x22)/
      (7.281820681114894e-58 + 6.793723008169782e-50*x + 1.0073238134263982e-42*x2 + 3.8567133904345664e-36*x3 + 4.6749055427591935e-30*x4 +
      1.9992235460663164e-24*x5 + 3.2233724058457452e-19*x6 + 2.051207320369606e-14*x7 + 5.334552198382988e-10*x8 + 5.833728132253219e-6*x9 +
      0.02748383972843651*x10 + 56.954350298361284*x11 + 52824.406590310646*x12 + 2.2217655338084057e7*x13 + 4.267549397728813e9*x14 +
      3.733806109621652e11*x15 + 1.459513002063948e13*x16 + 2.4180382199020853e14*x17 + 1.4786509784350255e15*x18 + 1.8740406611426415e15*x19 -
      3.345323820802959e15*x20 - 1.2075997985771218e15*x21 + 1.136213305508547e15*x22);

    if(!(aL >= 0.0001 && aL <= 1.0001))
    {
      aL = fmax(0.0001, fmin(aL, 1.0001));
    }

    aL2 = aL * aL;
    xi = 1. / aL2  -  1.;
  }

  double H200;
  double H240;

  if(fabs(aL - 1.) <= daL)  // taylor expansion around aL = 1
  {
    H200 = 0.0005624737680297483 + aL*(0.7797825941381623 + aL*(0.025694875450339216 +
         aL*(0.32086022134934217 + aL*(-0.19450407011549792 +
               aL*(0.09353199402114554 + aL*(-0.032113982358561236 + (0.006866621222216209 - 0.0006807274751769723*aL)*aL))))));


    H240 = -0.005087625772544818 + aL*(0.05695884534842821 + aL*(-0.3018911943816236 +
         aL*(1.019465776582752 + aL*(-0.5047550486743451 + aL*(0.16392338295963027 +
                  aL*(-0.030830985240682104 + (0.0020622455600826832 + 0.0001546036183062223*aL)*aL))))));
  }
  else if(aL > 0. && aL < 1. - daL) // 0 < aL < 1
  {
    double aL2 = aL * aL;
    double aL4 = aL2 * aL2;
    double sqrt_xi = sqrt(1./aL2  -  1.);
    double t = atan(sqrt_xi) / sqrt_xi;

    H200 = (aL2  +  t) / 2.;
    H240 = (2.  +  aL2  -  3. * t) * aL4 / (aL2 - 1.) / (aL2 - 1.);
  }
  else if(aL > 1. + daL)  // aL > 1
  {
    double aL2 = aL * aL;
    double aL4 = aL2 * aL2;
    double sqrt_xi = sqrt(1. - 1./aL2);
    double t = atanh(sqrt_xi) / sqrt_xi;

    H200 = (aL2  +  t) / 2.;
    H240 = (2.  +  aL2  -  3. * t) * aL4 / (aL2 - 1.) / (aL2 - 1.);
  }
  else
  {
    printf("compute_I240 error: aL = %lf is out of bounds", aL);
    exit(-1);
  }

  double I_240 = 0.5 * e * H240 / H200;

  return I_240;
}


double de_dt_ffe(double e, double pi, double t, double prefactor, double etas, double L04, double t0, double aL0, double z0, double w0)
{
  return (pi  -  4./3. * e) / t;
}


double compute_c0(double e, double pi, double t, double prefactor, double etas, double L04, double t0, double aL0)
{
  double p = e/3.;
  double T = pow(e / prefactor, 0.25);
  double taupi = 5. * etas / T;
  double x = t0 * aL0 / t;

  double e0 = 3. * L04 * HE_function(x) / (M_PI * M_PI);
  double pi0 = L04 * (HT_function(x)/2.  -  HL_function(x)) / (M_PI * M_PI);

  if(t == t0)
  {
    return 1.;
  }

  return pi * pi / (pi * pi0  + (e0 - e) * 16. * p / 15.);
}


double compute_c1(double e, double pi, double t, double prefactor, double etas, double L04, double t0, double aL0)
{
  double p = e/3.;
  double T = pow(e / prefactor, 0.25);
  double taupi = 5. * etas / T;
  double x = t0 * aL0 / t;

  double e0 = 3. * L04 * HE_function(x) / (M_PI * M_PI);
  double pi0 = L04 * (HT_function(x)/2.  -  HL_function(x)) / (M_PI * M_PI);

  if(t == t0)
  {
    return 0.;
  }

  return (t / taupi) * pi * (e0 - e) / (pi * pi0  + (e0 - e) * 16. * p / 15.);
}


double dpi_dt_ffe(double e, double pi, double t, double prefactor, double etas, double L04, double t0, double aL0, double z0, double w0)
{
  double p = e/3.;
  double T = pow(e / prefactor, 0.25);
  double taupi = 5. * etas / T;

  double x = t0 * aL0 / t;
  double A = compute_A(x, L04);

  double e0 = 3. * L04 * HE_function(x) / (M_PI * M_PI);
  double pi0 = L04 * (HT_function(x)/2.  -  HL_function(x)) / (M_PI * M_PI);

  double c0 = pi * pi / (pi * pi0  + (e0 - e) * 16. * p / 15.);
  double c1 = pi * (e0 - e) / (pi * pi0  + (e0 - e) * 16. * p / 15.);

  if(t == t0)
  {
    c0 = 1.;
    c1 = 0.;
  }

  if(std::isnan(c0) || std::isinf(c0))
  {
    printf("dpi_dt_ffe error: c0 = %lf\n", c0);
    //c0 = 1.;
  }

  if(std::isnan(c1) || std::isinf(c1))
  {
    printf("dpi_dt_ffe error: c1 = %lf\n", c1);
    //c1 = 0.;
  }

  double df0G = c0 * (A  -  31./15. * p);
  double df1G = - c1 * (608. * p  +  217. * pi) / 315.;

  return - pi / taupi  +  16. * p / (15. * t)  +  (df0G + df1G) / t;

}



void run_far_from_equilibrium_bjorken(double T0, double pl0, double t0, double dt, int Nt, double prefactor, double etas, double L0, double aL0)
{
  double t  = t0;
  double e = prefactor * pow(T0, 4);
  double pi = e/3.  -  pl0;
  double L04 = pow(L0, 4);
  double taupi0 = 5. * etas / T0;
  double w0 = t0 / taupi0;

  printf("Lambda_0 = %lf fm^-1\n", L0);
  printf("T_0      = %lf fm^-1\n", T0);
  printf("aL_0     = %lf\n\n", aL0);

  FILE *shear;
  // FILE *energy;
  // FILE *z0_file;
  // FILE *c0_file;
  // FILE *c1_file;

  shear = fopen("results/pibar_FFE.dat", "w");
  // energy = fopen("results/energy_FFE.dat", "w");
  // z0_file = fopen("results/z0_FFE.dat", "w");
  // c0_file = fopen("results/c0_FFE.dat", "w");
  // c1_file = fopen("results/c1_FFE.dat", "w");

  double z0 = 0;

  for(int i = 1; i < Nt; i++)
  {
    double TL = pow(e / prefactor, 0.25);
    double taupi_inverse_L = TL / (5. * etas);                       // left point

    double c0 = compute_c0(e, pi, t, prefactor, etas, L04, t0, aL0);
    double c1 = compute_c1(e, pi, t, prefactor, etas, L04, t0, aL0);

    fprintf(shear, "%.6e\t%.6e\n", t, 3. * pi / e);                  // pi / peq
    // fprintf(energy, "%.6e\t%.6e\n", t, e);                           // e
    // fprintf(z0_file, "%.6e\t%.6e\n", t, z0);                         // z0 integration via trapezoid rule
    // fprintf(c0_file, "%.6e\t%.6e\n", t, c0);                         // c0 coefficient in FFE
    // fprintf(c1_file, "%.6e\t%.6e\n", t, c1);                         // c1 coefficient in FFE

    double de1  = dt *  de_dt_ffe(e, pi, t, prefactor, etas, L04, t0, aL0, z0, w0);
    double dpi1 = dt * dpi_dt_ffe(e, pi, t, prefactor, etas, L04, t0, aL0, z0, w0);

    double TR = pow((e + de1) / prefactor, 0.25);                   // intermediate right point
    double taupi_inverse_R = TR / (5. * etas);

    double dz01 = dt * (taupi_inverse_L + taupi_inverse_R) / 2.;    // intermediate iteration

    double de2  = dt *  de_dt_ffe(e + de1, pi + dpi1, t + dt, prefactor, etas, L04, t0, aL0, z0 + dz01, w0);
    double dpi2 = dt * dpi_dt_ffe(e + de1, pi + dpi1, t + dt, prefactor, etas, L04, t0, aL0, z0 + dz01, w0);

    e  += (de1  + de2)  / 2.;                                       // RK2
    pi += (dpi1 + dpi2) / 2.;

    TR = pow(e / prefactor, 0.25);                                  // recompute right point
    taupi_inverse_R = TR / (5. * etas);

    z0 += dt * (taupi_inverse_L + taupi_inverse_R) / 2.;            // update z0 via trapezoid rule integration

    t += dt;
  }

  fclose(shear);
  // fclose(energy);
  // fclose(z0_file);
  // fclose(c0_file);
  // fclose(c1_file);
}

