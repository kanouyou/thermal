#ifndef IFdmNbTi_HH
#define IFdmNbTi_HH

namespace FDM
{ class IFdmNbTi; }

class FDM::IFdmNbTi
{
  public:
    static IFdmNbTi* Instance();

    ~IFdmNbTi();

    void SetField(const double B) { fField = B; }

    double GetCapacity(double T);

    double GetDensity() const { return 6538.; }

    double GetCriticalCurrent(double T);

  private:
    static IFdmNbTi* inst;

    IFdmNbTi();

  private:
    double fField;
};

#endif
