
#ifndef VISCOUSBJORKEN_H_
#define VISCOUSBJORKEN_H_

void run_second_order_viscous_bjorken(double T0, double pl0, double t0, double dt, int Nt, double prefactor, double etas);
void run_third_order_viscous_bjorken(double T0, double pl0, double t0, double dt, int Nt, double prefactor, double etas);

#endif