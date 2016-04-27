#include <iostream>
#include <cmath>
#include "IFdmKapton.h"
#include "IFdmException.h"

FDM::IFdmKapton* FDM::IFdmKapton::inst = NULL;

FDM::IFdmKapton::IFdmKapton()
    : fT0(4.2), fTf(273.)
{}

FDM::IFdmKapton::~IFdmKapton()
{}

FDM::IFdmKapton* FDM::IFdmKapton::Instance()
{
  if (inst==NULL)
    inst = new FDM::IFdmKapton();

  return inst;
}

void FDM::IFdmKapton::SetTemperatureRange(double T0, double Tf)
{
  if (T0>Tf) {
    FDM::IFdmException except("Tmin is larger than Tmax!", "IFdmKapton", "SetTemperatureRange()");
    throw except;
  }

  fT0 = T0;
  fTf = Tf;
}

double FDM::IFdmKapton::GetThermalConductivity(double T)
{
  // fitting parameter
  const int n = 8;
  const double p[n] = {5.73101, -39.5199, 79.9313, -83.8572, 50.9157,
                       -17.9835, 3.42413, -0.27133};

  //
  double k;
  double ax = 0.;

  if (T<4.3)
    k = 0.0378 + 0.00161*T;
  else {
    for (int i=0; i<n; i++)
      ax += p[i] * pow(log10(T),i);
    //
    k = pow(10, ax);
  }

  return k;
}

double FDM::IFdmKapton::GetCapacity(double T)
{
  // fitting parameters
  const int n = 8;
  const double p[n] = {-1.3684,  0.65892, 2.8719, 0.42651, -3.0088,
                        1.9558, -0.51998, 0.051574};

  //
  double ax = 0.;
  for (int i=0; i<n; i++)
    ax += p[i] * pow(log10(T),i);

  double C = pow(10, ax);

  return C;
}

TGraph* FDM::IFdmKapton::DrawCapacity()
{
  TGraph* gr = new TGraph();

  const int    nT = 500;
  const double dT = (fTf-fT0) / nT;

  //
  double T, C;

  for (int i=0; i<nT; i++) {
    T = fT0 + i * dT;
    C = GetCapacity(T);
    gr->SetPoint(i, T, C);
  }

  gr->SetLineColor(kSpring+9);
  gr->SetLineWidth(2);
  gr->SetTitle("Kapton; Temperature [K]; Heat Capacity [W/m/K]");

  return gr;
}
