#ifndef IFdmCopper_HH
#define IFdmCopper_HH

#include <TGraph.h>

namespace FDM
{ class IFdmCopper; }

class FDM::IFdmCopper
{
  public:
    static IFdmCopper* Instance();

    ~IFdmCopper();

    void SetRRR(const double RRR) { fRRR = RRR; }

    void SetField(double B);

    void SetTemperatureRange(double T0, double Tf);

    double GetResistivity(double T);

    double GetThermalConductivity(double T);

    double GetCapacity(double T);

    double GetDensity() const { return 8960.; }

    TGraph* DrawCapacity();

    TGraph* DrawResistivity();

  private:
    static IFdmCopper* inst;

    IFdmCopper();

  private:
    double fRRR;
    double fField;
    double fT0;
    double fTf;
};

#endif
