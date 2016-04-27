#ifndef IFdmHeatContainer_HH
#define IFdmHeatContainer_HH

namespace FDM
{ class IFdmHeatContainer; }

class FDM::IFdmHeatContainer
{
  public:
    IFdmHeatContainer();

    ~IFdmHeatContainer();

    void SetId(int x, int y, int z);

    void SetPosition(int x, int y, int z);

    void SetTime(double t) { fTime = t; }

    void SetPreTemperature(double T) { fPreT = T; }

    void SetPostTemperature(double T) { fPostT = T; }

    void SetCapacity(double C) { fCapacity = C; }

    void SetConductivity(double kx, double ky, double kz);

    void SetDensity(double rho) { fDensity = rho; }

    void SetLength(double dx, double dy, double dz);

    void SetHeat(double qx, double qy, double qz);

    void SetHeatSource(double Q) { fSource = Q; }

    void SetRRR(double RRR) { fRRR = RRR; }

    void SetField(double B) { fField = B; }

    int* GetId() { return fId; }

    double* GetPosition() { return fPos; }

    double GetTime() { return fTime; }

    double GetPreTemperature() { return fPreT; }

    double GetPostTemperature() { return fPostT; }

    double GetCapacity() { return fCapacity; }

    double* GetConductivity() { return fConductivity; }

    double GetDensity() { return fDensity; }

    double* GetLength() { return fLength; }

    double* GetHeat() { return fQ; }

    double GetHeatSource() { return fSource; }

    double GetRRR() { return fRRR; }

    double GetField() { return fField; }


  private:
    int* fId;
    double* fPos;
    double fPreT;
    double fPostT;
    double fCapacity;
    double* fConductivity;
    double fDensity;
    double* fLength;
    double* fQ; 
    double  fSource;
    double  fTime;
    double  fRRR;
    double  fField;
};

#endif
