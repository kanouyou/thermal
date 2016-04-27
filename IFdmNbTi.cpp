#include <iostream>
#include <cmath>
#include "IFdmNbTi.h"

FDM::IFdmNbTi* FDM::IFdmNbTi::inst = NULL;

FDM::IFdmNbTi::IFdmNbTi()
    : fField(0.)
{}

FDM::IFdmNbTi::~IFdmNbTi()
{}

FDM::IFdmNbTi* FDM::IFdmNbTi::Instance()
{
  if (inst==NULL)
    inst = new FDM::IFdmNbTi();

  return inst;
}

double FDM::IFdmNbTi::GetCapacity(double T)
{
  // fitting parameter
  const int n = 5;
  const double p1[n] = {0.0, 64.*fField, 0.0, 49.1, 0.0};
  const double p2[n] = {0.0, 928., 0.0, 16.24, 0.0};
  const double p3[n] = {41383., -7846.1, 553.71, 11.9838, -0.2177};
  const double p4[n] = {-1.53e+6, 83022., -716.3, 2.976, -0.00482};
  const double p5[n] = {1.24e+6, 13706., -51.66, 0.09296, -6.29e-5};
  const double p6[n] = {2.45e+6, 955.5, -0.257, 0., 0.};

  //
  const double Tc  = 9.4;
  const double rho = 6538.;

  //
  double C = 0.;

  if (T>0. && T<Tc) {
    for (int i=0; i<n; i++)
      C += p1[i] * pow(T, i);    
  }
  else if (T>=Tc && T<28.358) {
    for (int i=0; i<n; i++)
      C += p2[i] * pow(T, i);
  }
  else if (T>=28.358 && T<50.99) {
    for (int i=0; i<n; i++)
      C += p3[i] * pow(T, i);
  }
  else if (T>=50.99 && T<165.8) {
    for (int i=0; i<n; i++)
      C += p4[i] * pow(T, i);
  }
  else if (T>=165.8 && T<496.54) {
    for (int i=0; i<n; i++)
      C += p5[i] * pow(T, i);
  }
  else if (T>=496.54) {
    for (int i=0; i<n; i++)
      C += p6[i] * pow(T, i);
  }

  return C / rho;
}

double FDM::IFdmNbTi::GetCriticalCurrent(double T)
{
  // fitting equation from L. Bottura's paper
  const int ni = 5;
  const double n = 1.7;
  const double Tc0 [ni] = {9.2, 8.5, 8.9, 9.2, 9.35};
  const double Bc20[ni] = {14.5, 14.2, 14.4, 14.4, 14.25};
  const double C0  [ni] = {23.8, 28.6, 28.5, 37.7, 28.4};
  const double alp [ni] = {0.57, 0.76, 0.64, 0.89, 0.80};
  const double beta[ni] = {0.90, 0.85, 0.75, 1.10, 0.89};
  const double gam [ni] = {1.90, 1.76, 2.30, 2.09, 1.87};

  //
  const int m = 0;
  double t   = T / Tc0[m];
  double Bc2 = Bc20[m] * (1 - pow(t,n));
  double b   = fField / Bc2;

  // normalized critical current density
  double Jc = C0[m] * pow(b,alp[m]) * pow((1-b),beta[m]) * pow((1-pow(t,n)),gam[m]) / fField;
  const double J0 = 3000.;    // [A/mm2]

  Jc = J0 * Jc;
  return Jc;
}
