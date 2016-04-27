#include <iostream>
#include <cmath>
#include "IFdmAluminium.h"
#include "IFdmException.h"

FDM::IFdmAluminium* FDM::IFdmAluminium::inst = NULL;

FDM::IFdmAluminium* FDM::IFdmAluminium::Instance()
{
  if (inst==NULL)
    inst = new FDM::IFdmAluminium();

  return inst;
}

FDM::IFdmAluminium::IFdmAluminium()
    : fRRR(1.), fField(1.),
      fT0(4.2), fTf(273.)
{}

FDM::IFdmAluminium::~IFdmAluminium()
{}

void FDM::IFdmAluminium::SetRRR(const double RRR)
{
  fRRR = RRR;
}

void FDM::IFdmAluminium::SetField(const double B)
{
  fField = B;
}

void FDM::IFdmAluminium::SetTemperatureRange(double T0, double Tf)
{
  if (Tf<T0) {
    FDM::IFdmException except("Tmax is less than Tmin!", "IFdmAluminium", "SetTemperatureRange()");
    throw except;
  }

  fT0 = T0;
  fTf = Tf;
}

double FDM::IFdmAluminium::GetResistivity(double T)
{
  double rho0, rhoi, rhoi0, rhoB;

  // fitting parameter
  const int n = 7;
  const double p[n] = {1.671e-17, 4.36, 2.841e+10, 1.18, 64.0,
                       4.428, 1.2031};

  // fitting parameter for Al magnetoresistance
  const int nb = 5;
  const double pb[nb] = {3.62857, 2.90419e-5, 3.79649e+6, 10975.9, 0.761609};

  // resistivity at room temperature [Ohm*m]
  const double rhoRT = 2.75e-8;

  rho0  = rhoRT / fRRR;
  rhoi  = p[0] * pow(T,p[1]) / (1 + p[0]*p[2]*pow(T,p[1]-p[3])*exp(-pow(p[4]/T,p[5])));
  rhoi0 = p[6] * rhoi * rho0 / (rho0 + rhoi);

  // resistivity
  double rho = rho0 + rhoi + rhoi0;

  // magnetoresistance
  if (fField>0.) {
    double h = fField * 10. * rhoRT / rho;
    rhoB = h*h*(pb[0] - pb[1]*h)*rho / (pb[2] + pb[3]*h + pb[4]*h*h) + rho;
    rho  = rhoB;
  }

  return rho;
}

double FDM::IFdmAluminium::GetThermalConductivity(double T)
{
  const double Lwf = 2.44e-8;
  double rho = GetResistivity(T);
  double k = Lwf * T / rho;

  return k;
}

double FDM::IFdmAluminium::GetCapacity(double T)
{
  // switch point
  double sp1 = 22.67;
  double sp2 = 46.0;

  // fitting parameters
  const int n = 4;
  const double p1[n] = {-0.207489, 0.165759, -0.0142572, 0.00146459};
  const double p2[n] = {7.88e-13, 6.93201, -0.07139, 46.4363};
  const double p3[n] = {6.273517, -0.5469, 0.000925, -156.932};

  //
  double C;
  if (T>0. && T<sp1)
    C = p1[0] + p1[1]*T + p1[2]*pow(T,2) + p1[3]*pow(T,3);
  else if (T>=sp1 && T<sp2)
    C = p2[0] * pow(T,p2[1]) * exp(p2[2]*T) * exp(p2[3]/T) * 4.186e+3;
  else if (T>=sp2)
    C = p3[0] * pow(T,p3[1]) * exp(p3[2]*T) * exp(p3[3]/T) * 4.186e+3;
  else {
    FDM::IFdmException except("Temperature is out of range!", "IFdmAluminium", "GetCapacity()");
    throw except;
  }

  return C;
}

TGraph* FDM::IFdmAluminium::DrawCapacity()
{
  TGraph* gr = new TGraph();

  const int nT = 500;
  double dT = (fTf - fT0) / nT;
  double T, C;

  for (int i=0; i<nT; i++) {
    T = fT0 + dT * i;
    C = GetCapacity(T);
    gr->SetPoint(i, T, C);
  }

  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->SetTitle("Aluminium; Temperature [K]; Heat Capacity [J/kg/K]");

  return gr;
}

TGraph* FDM::IFdmAluminium::DrawResistivity()
{
  TGraph* gr = new TGraph();

  const int nT = 500;
  double dT = (fTf - fT0) / nT;
  double T, rho;

  for (int i=0; i<nT; i++) {
    T   = fT0 + dT * i;
    rho = GetResistivity(T);
    gr->SetPoint(i, T, rho);
  }

  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->SetTitle("Aluminium; Temperature [K]; Electric Resistivity [#Omega #cdot m]");

  return gr;
}
