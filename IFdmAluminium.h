#ifndef IFdmAluminium_HH
#define IFdmAluminium_HH

#include <TGraph.h>

namespace FDM
{ class IFdmAluminium; }

class FDM::IFdmAluminium
{
  public:
    static IFdmAluminium* Instance();

    ~IFdmAluminium();

    void SetRRR(const double RRR);

    void SetField(const double B);

    void SetTemperatureRange(double T0, double Tf);

    double GetDensity() { return 2700.; }

    double GetResistivity(double T);

    double GetThermalConductivity(double T);

    double GetCapacity(double T);

    TGraph* DrawCapacity();

    TGraph* DrawResistivity();

  private:
    static IFdmAluminium* inst;

    IFdmAluminium();

  private:
    double fRRR;
    double fField;
    double fT0;
    double fTf;
};

#endif
