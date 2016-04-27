#ifndef IFdmKapton_HH
#define IFdmKapton_HH

#include <TGraph.h>

namespace FDM
{ class IFdmKapton; }

class FDM::IFdmKapton
{
  public:
    static FDM::IFdmKapton* Instance();

    ~IFdmKapton();

    void SetTemperatureRange(double T0, double Tf);

    double GetThermalConductivity(double T);

    double GetCapacity(double T);

    double GetDensity() const { return 1420.; }

    TGraph* DrawCapacity();

  private:
    static FDM::IFdmKapton* inst;

    IFdmKapton();

  private:
    double fT0;
    double fTf;

};

#endif
