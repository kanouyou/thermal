#include <iostream>
#include <fstream>
#include "IFdmException.h"
#include "IFdmRadiation.h"

using std::string;
using std::vector;

FDM::IFdmRadiation::IFdmRadiation()
    : fMesh(NULL), cnt(0), fTime(0.)
{}

FDM::IFdmRadiation::~IFdmRadiation()
{
  if (fMesh) delete[] fMesh;
}

void FDM::IFdmRadiation::Load(const string& file)
{
  std::ifstream input(file.c_str()); 

  int    ibuff[3];
  double dbuff[2];
  IRadiationMap* radmap = NULL;

  if (!input) {
    FDM::IFdmException except("File is not founded!", "IFdmRadiation", "Load");
    throw except;
  }
  else {
    while (true) {
      input >> ibuff[0] >> ibuff[1] >> ibuff[2] >> dbuff[0] >> dbuff[1];
      if (!input) break;
      radmap = new IRadiationMap();
      radmap->SetXYZ(ibuff[0], ibuff[1], ibuff[2]);
      radmap->SetFlux(dbuff[0]);
      radmap->SetHeat(dbuff[1]);
      fData[cnt] = radmap;
      cnt++;
    }
  }
}

void FDM::IFdmRadiation::SetMesh(const int mz, const int mp, const int mr)
{
  if (!fMesh) fMesh = new int[3];

  fMesh[0] = mz;
  fMesh[1] = mp;
  fMesh[2] = mr;
}

double FDM::IFdmRadiation::GetRRR(int ii, int jj, int kk)
{
  ii--; jj--; kk--;

  const int dnz = fData[cnt-1]->X() + 1;
  const int dnp = fData[cnt-1]->Y() + 1;
  const int dnr = fData[cnt-1]->Z() + 1;

  const int dz = fMesh[0] / dnz;
  const int dp = fMesh[1] / dnp;
  const int dr = fMesh[2] / dnr;

  const int iz = ii / dz;
  const int ip = jj / dp;
  const int ir = kk / dr;

  double flux = 0.;

  for (std::map<int, IRadiationMap*>::iterator it=fData.begin(); it!=fData.end(); ++it) {
    if (it->second->X()==iz && it->second->Y()==ip && dnr-it->second->Z()-1==ir)
      flux = it->second->GetFlux();
  }

  double timefactor = fTime / 365.;
  double neufactor  = 2.7e-22;
  double RRR = 0.;

  if (kk%2==0)
    // conductor
    RRR = 27. / (0.0135 + flux * timefactor * neufactor);
  else if (kk%2==1 && kk<19)
    // strip
    RRR = 27. / (0.0675 + flux * timefactor * neufactor);

  return RRR;
}

double FDM::IFdmRadiation::GetHeat(int ii, int jj, int kk)
{
  ii--; jj--; kk--;

  const int dnz = fData[cnt-1]->X() + 1;
  const int dnp = fData[cnt-1]->Y() + 1;
  const int dnr = fData[cnt-1]->Z() + 1;

  const int dz = fMesh[0] / dnz;
  const int dp = fMesh[1] / dnp;
  const int dr = fMesh[2] / dnr;

  const int iz = ii / dz;
  const int ip = jj / dp;
  const int ir = kk / dr;

  double heat = 0.;

  for (std::map<int, IRadiationMap*>::iterator it=fData.begin(); it!=fData.end(); ++it) {
    if (it->second->X()==iz && it->second->Y()==ip && dnr-it->second->Z()-1==ir)
      heat = it->second->GetHeat();
  }

  return heat;
}


