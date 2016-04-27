#include <iostream>
#include <cmath>
#include "IFdmCopper.h"
#include "IFdmException.h"

FDM::IFdmCopper* FDM::IFdmCopper::inst = NULL;

FDM::IFdmCopper::IFdmCopper()
    : fRRR(1.), fField(1.),
      fT0(4.2), fTf(273.)
{}

FDM::IFdmCopper::~IFdmCopper()
{}

FDM::IFdmCopper* FDM::IFdmCopper::Instance()
{
  if (inst==NULL)
    inst = new FDM::IFdmCopper();

  return inst;
}

void FDM::IFdmCopper::SetField(double B)
{
  if (B<0)
    B *= -1.;

  std::cerr << "Note: your inputed field is negative!" << std::endl;

  fField = B;
}

void FDM::IFdmCopper::SetTemperatureRange(double T0, double Tf)
{
  if (Tf < T0) {
    FDM::IFdmException except("Tmax is less than Tmin!", "IFdmCopper", "SetTemperatureRange()");
    throw except;
  }

  fT0 = T0;
  fTf = Tf;
}

double FDM::IFdmCopper::GetResistivity(double T) 
{
  // fitting parameter
  const int n = 7;
  const double p[n] = {1.171e-17, 4.49, 3.841e+10, 1.14, 50.0,
                       6.428, 0.4531};

  // resistivity at room temperature [Ohm*m]
  double rhoRT = 1.553e-8;
  double rho0  = rhoRT / fRRR;
  double rhoi  = p[0] * pow(T,p[1]) / (1 + p[0] * p[2] * pow(T, p[1]-p[3]) * exp(-pow(p[4]/T,p[5])));
  double rhoi0 = p[6] * rhoi * rho0 / (rho0 + rhoi);

  double rho   = rho0 + rhoi + rhoi0;

  // megnetoresistance
  if (fField>0.) {
    // fitting parameter
    const int nb = 5;
    const double pb[nb] = {-2.662, 0.3168, 0.6229, -0.1839, 0.01827};
    //
    double ax = pb[0]*pow(log10(fRRR*fField),0) + pb[1]*pow(log10(fRRR*fField),1) 
              + pb[2]*pow(log10(fRRR*fField),2) + pb[3]*pow(log10(fRRR*fField),3) 
              + pb[4]*pow(log10(fRRR*fField),4); 
    double rhoB = rho * (1 + pow(10,ax));
    rho = rhoB;
  }

  return rho;
}

double FDM::IFdmCopper::GetThermalConductivity(double T)
{
  double beta  = 0.634 / fRRR;
  double betar = beta / 0.0003;

  // fitting parameter
  const int n = 7;
  const double p[n] = {1.754e-8, 2.763, 1102.0, -0.165, 70.0,
                       1.756, 0.838/pow(betar,0.1661)};

  //
  double W0  = beta / T;
  double Wi  = p[0]*pow(T,p[1]) / (1 + p[0]*p[2]*pow(T,(p[1]+p[3]))*exp(-pow(p[4]/T,p[5])));
  double Wi0 = p[6] * Wi * W0 / (Wi + W0);
  double k   = 1. / (W0 + Wi + Wi0);

  const double B = fField;

  // magnetoconductivity
  if (fField>0.) {
    double rho = GetResistivity(T);
    SetField(0.);
    k = GetResistivity(T) * k / rho;
  }

  SetField(B);

  return k;
};

double FDM::IFdmCopper::GetCapacity(double T)
{
  // fitting parameter
  const int n = 8;
  const double p[n] = {-1.91844, -0.15973, 8.61013, -18.996,
                        21.9661, -12.7328, 3.54322, -0.3797};

  //
  double ax = 0.;
  for (int i=0; i<n; i++)
    ax += p[i] * pow(log10(T), i);

  double C = pow(10, ax);

  return C;
}

TGraph* FDM::IFdmCopper::DrawCapacity()
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

  gr->SetLineColor(kAzure+4);
  gr->SetLineWidth(2);
  gr->SetTitle("Copper; Temperature [K]; Heat Capacity [J/kg/K]");

  return gr;
}

TGraph* FDM::IFdmCopper::DrawResistivity()
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

  gr->SetLineColor(kAzure+4);
  gr->SetLineWidth(2);
  gr->SetTitle("Copper; Temperature [K]; Electric Resistivity [#Omega #cdot m]");

  return gr;
}
