#ifndef IFdmRadiation_HH
#define IFdmRadiation_HH

#include <string>
#include <vector>
#include <map>

// *******************************
namespace FDM
{ class IRadiationMap; }

class FDM::IRadiationMap
{
  public:
    IRadiationMap() {}

    ~IRadiationMap() {}

    void SetXYZ(int z, int p, int r) {
      fZ = z;
      fP = p;
      fR = r;
    }

    void SetFlux(double flux) { fNeu = flux; }

    void SetHeat(double heat) { fHeat = heat; }

    int X() { return fZ; }

    int Y() { return fP; }

    int Z() { return fR; }

    double GetFlux() { return fNeu; }

    double GetHeat() { return fHeat; }

  private:
    int    fZ;
    int    fP;
    int    fR;
    double fNeu;
    double fHeat;
};



namespace FDM
{ class IFdmRadiation; }

class FDM::IFdmRadiation
{
  public:
    IFdmRadiation();

    ~IFdmRadiation();

    void Load(const std::string& file);

    void SetMesh(const int mz, const int mp, const int mr);

    void SetOperTime(double time) { fTime = time; }

    double GetRRR(int, int, int);

    double GetHeat(int, int, int);

  private:
    int* fMesh;
    int  cnt;
    double fTime;
    std::map<int, FDM::IRadiationMap*> fData;
};


#endif
