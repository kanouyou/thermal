#include <iostream>
#include "IFdmHeatContainer.h"

FDM::IFdmHeatContainer::IFdmHeatContainer()
    : fId(NULL), fPos(NULL),
      fPreT(0.), fPostT(0.),
      fCapacity(0.), fConductivity(NULL),
      fDensity(0.), fLength(NULL),
      fQ(NULL), fSource(0.),
      fTime(0.), fRRR(0.),
      fField(0.)
{
  fConductivity = new double[3];
  fConductivity[0] = 0.;
  fConductivity[1] = 0.;
  fConductivity[2] = 0.;
}

FDM::IFdmHeatContainer::~IFdmHeatContainer()
{
  delete[] fId;
  delete[] fPos;
  delete[] fConductivity;
  delete[] fLength;
  delete[] fQ;
}

void FDM::IFdmHeatContainer::SetId(int x, int y, int z)
{
  if (fId==NULL)
    fId = new int[3];

  fId[0] = x;
  fId[1] = y;
  fId[2] = z;
}

void FDM::IFdmHeatContainer::SetPosition(int x, int y, int z)
{
  if (fPos==NULL)
    fPos = new double[3];

  fPos[0] = x;
  fPos[1] = y;
  fPos[2] = z;
}

void FDM::IFdmHeatContainer::SetConductivity(double kx, double ky, double kz)
{
  if (fConductivity==NULL)
    fConductivity = new double[3];

  fConductivity[0] = kx;
  fConductivity[1] = ky;
  fConductivity[2] = kz;
}

void FDM::IFdmHeatContainer::SetLength(double dx, double dy, double dz)
{
  if (fLength==NULL)
    fLength = new double[3];

  fLength[0] = dx;
  fLength[1] = dy;
  fLength[2] = dz;
}

void FDM::IFdmHeatContainer::SetHeat(double qx, double qy, double qz)
{
  if (fQ==NULL)
    fQ = new double[3];

  fQ[0] = qx;
  fQ[1] = qy;
  fQ[2] = qz;
}

