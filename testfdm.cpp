#include <iostream>
#include <string>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include "IFdmRadiation.h"
#include "IFdmAluminium.h"
#include "IFdmCopper.h"
#include "IFdmKapton.h"
#include "IFdmUnits.h"
#include "IFdmException.h"
#include "IFdmHeatContainer.h"

using namespace std;
using namespace FDM;

double GetAvgkz(double l1, double k1, double l2, double k2)
{
  double k = (l1 + l2) / (l1/k1 + l2/k2);
  return k;
}

double GetAvgkp(double A1, double k1, double A2, double k2)
{
  double k = (A1*k1 + A2*k2) / (A1 + A2);
  return k;
}

int main(int argc, char** argv)
{

  const char* filename = "./output/result1.txt";
  const char* matfile  = "./output/material1.txt";

  // mesh
  const int Mz = 270 + 2;
  const int Mp = 4 + 2;
  const int Mr = 19 + 2;

  IFdmRadiation* rad = new IFdmRadiation();
  rad->Load("./data/phitshield.txt");
  rad->SetOperTime(20.);
  rad->SetMesh(Mz, Mp, Mr);

  IFdmAluminium* al  = IFdmAluminium::Instance();
  IFdmKapton*    kap = IFdmKapton::Instance();

  al->SetRRR(200.);
  al->SetField(5.5);

  // aluminium strip
  const double lstrip = 1.*mm;
  const double ledge  = 10.*cm;
  const double rhoAl  = 2700.;
  const double CAl    = al->GetCapacity(4.5);
  double kAl  = al->GetThermalConductivity(4.5);
  // conductor
  const double lcdt   = 15.*mm;
  const double wcdt   = 4.73*mm;
  const double rhoCdt = 4000.;
  const double CCdt   = 0.16;
  al->SetRRR(80.);
  double kCdt = al->GetThermalConductivity(4.5);
  // insulation
  const double ltape  = 0.4*mm;
  const double wtape  = 0.8*mm;
  //const double rhoTape= 1420.;
  //const double CTape  = kap->GetCapacity(4.5);
  const double kTape  = kap->GetThermalConductivity(4.5);
  // resin
  const double wresin = 2*mm;
  const double kresin = 0.05;
  // shell
  const double lshell = 50.*mm;
  al->SetRRR(5.);
  const double kshell = al->GetThermalConductivity(4.5);

  //const double lz = (ltape*2 + wcdt) * 270;
  const double lp = 2.*M_PI*800*mm;
  //const double lr = (lcdt + 2*ltape + lstrip)*9 + lshell;

  //const double dz = lz / Mz;
  const double dp = lp / Mp;
  //const double dr = lr / Mr;

  const double t0 = 0.*sec;
  const double tf = 40*sec;
  const double dt = 0.004*msec;
  const int    nt = (tf - t0) / dt;

  // cell parameters
  const double hz  = ltape*2 + wcdt;
  const double hp  = dp;
  const double hr  = (ltape*2 + lcdt)/2 + lstrip/2;

  //const double kz  = GetAvgkz(wcdt, kCdt, wtape, kTape);
  //const double kp  = GetAvgkp(lcdt*wcdt, kCdt, (lcdt+ltape)*(wcdt+ltape)-lcdt*wcdt, kTape);
  //const double kr  = GetAvgkz(lcdt, kCdt, wtape, kTape);

  double kz, kp, kr;
  double T0 = 4.5;
  double B0 = 5.5;

  ofstream output;
  output.open(filename, ios::out);

  ofstream outmat;
  outmat.open(matfile, ios::out);

  // initialization
  static IFdmHeatContainer* hc[Mz*Mp*Mr];
  int ic;

  for (int i=0; i<Mz; i++) {
    for (int j=0; j<Mp; j++) {
      for (int k=0; k<Mr; k++) {
        ic = i*Mp*Mr + j*Mr + k;
        hc[ic] = new IFdmHeatContainer();
        hc[ic]->SetId(i, j, k);
        hc[ic]->SetPreTemperature(4.5);
        hc[ic]->SetPostTemperature(4.5);
        hc[ic]->SetHeatSource(rad->GetHeat(i,j,k));

        al->SetRRR(rad->GetRRR(i,j,k));
        al->SetField(B0);
        hc[ic]->SetRRR(rad->GetRRR(i,j,k));
        hc[ic]->SetField(B0);

        kCdt = al->GetThermalConductivity(T0);
        kAl  = al->GetThermalConductivity(T0);

        kz = GetAvgkz(wcdt, kCdt, wtape, kTape);
        kp = GetAvgkp(lcdt*wcdt, kCdt, (lcdt+ltape)*(wcdt+ltape)-lcdt*wcdt, kTape);
        kr = GetAvgkz(lcdt, kCdt, wtape, kTape);

        switch (k) {
          // for boundary
          case 0:
            hc[ic]->SetLength(lstrip/2, hp, lstrip/2);
          // aluminium strip
          case 1:
          case 3:
          case 5:
          case 7:
          case 9:
          case 11:
          case 13:
          case 15:
          case 17:
            hc[ic]->SetConductivity(kAl, kresin, kAl);
            hc[ic]->SetDensity(rhoAl);
            hc[ic]->SetCapacity(CAl);
            hc[ic]->SetLength(hz, hp, hr);
            if (i==0 || i>=Mz-2)
              hc[ic]->SetLength(hz/2+ledge/2, hp, hr);
            break;
          // conductor
          case 2:
          case 4:
          case 6:
          case 8:
          case 10:
          case 12:
          case 14:
          case 16:
          case 18:
            hc[ic]->SetConductivity(kz, kp, kr);
            hc[ic]->SetDensity(rhoCdt);
            hc[ic]->SetCapacity(CCdt);
            hc[ic]->SetLength(hz, hp, hr);
            if (i==0 || i>=Mz-2)
              hc[ic]->SetLength(hz/2, hp, (ltape*2 + lcdt)/2 + lshell/2);
            if (k==18)
              hc[ic]->SetLength(hz, hp, (ltape*2 + lcdt)/2 + lshell/2);
            break;
          // shell
          case 19:
            hc[ic]->SetConductivity(kshell, kshell, GetAvgkz(wresin,kresin,lshell,kshell));
            hc[ic]->SetDensity(rhoAl);
            hc[ic]->SetCapacity(CAl);
            hc[ic]->SetLength(hz, hp, lshell/2);
            hc[ic]->SetRRR(5.);
            if (i==0 || i>=Mz-2)
              hc[ic]->SetLength(hz/2, hp, lshell/2);
            break;
          default:
            break;
        }
      }
    }
  }
  
  double Qxyz[Mz*Mp*Mr];
  double QQ  [Mz*Mp*Mr];
  float  T   [Mz*Mp*Mr];
  float  preT[Mz*Mp*Mr];
  float  time; 
  double kxyz[3]; 
  double hi  [3]; 
  double hi_1[3]; 
  double RadQ, C, rho, Q;


  for (int i=0; i<Mz; i++) {
    for (int j=0; j<Mp; j++) {
      for (int k=0; k<Mr; k++) {
        preT[i*Mp*Mr + j*Mr + k] = 4.5;
        T   [i*Mp*Mr + j*Mr + k] = 4.5;
      }
    }
  }

  for (int i=1; i<Mz-1; i++) {
    for (int j=1; j<Mp-1; j++) {
      for (int k=1; k<Mr-1; k++) {
        ic = i*Mp*Mr + j*Mr + k;
        // *************************************************
        // MATERIAL PROPERTY:
        outmat << setw(5) << left << i
               << setw(5) << left << j
               << setw(5) << left << k
               << setw(17) << right << scientific << hc[ic]->GetDensity()
               << setw(17) << right << scientific << hc[ic]->GetCapacity()
               << setw(17) << right << scientific << hc[ic]->GetRRR()
               << setw(17) << right << scientific << hc[ic]->GetConductivity()[0]
               << setw(17) << right << scientific << hc[ic]->GetConductivity()[1]
               << setw(17) << right << scientific << hc[ic]->GetConductivity()[2]
               << setw(17) << right << scientific << hc[ic]->GetHeatSource() << "\n";
      }
    }
  }
  outmat.close();


  for (int it=0; it<nt; it++) {
    time = t0 + dt*it;

    for (int i=0; i<Mz; i++) {
      for (int j=0; j<Mp; j++) {
        for (int k=0; k<Mr; k++) {
          preT[i*Mp*Mr + j*Mr + k] = T[i*Mp*Mr + j*Mr + k];
        }
      }
    }

    for (int i=1; i<Mz-1; i++) {
      for (int j=1; j<Mp-1; j++) {
        for (int k=1; k<Mr-1; k++) {
          kxyz[0] = hc[i*Mp*Mr + j*Mr + k]->GetConductivity()[0];
          kxyz[1] = hc[i*Mp*Mr + j*Mr + k]->GetConductivity()[1];
          kxyz[2] = hc[i*Mp*Mr + j*Mr + k]->GetConductivity()[2];
          hi  [0] = hc[i*Mp*Mr + j*Mr + k]->GetLength()[0];
          hi  [1] = hc[i*Mp*Mr + j*Mr + k]->GetLength()[1];
          hi  [2] = hc[i*Mp*Mr + j*Mr + k]->GetLength()[2];
          hi_1[0] = hc[(i-1)*Mp*Mr + j*Mr + k]->GetLength()[0];
          hi_1[1] = hc[i*Mp*Mr + (j-1)*Mr + k]->GetLength()[1];
          hi_1[2] = hc[i*Mp*Mr + j*Mr + (k-1)]->GetLength()[2];

/*
          Qxyz[0] = kxyz[0] * ((preT[i*Mp*Mr + j*Mr + k] - preT[(i-1)*Mp*Mr + j*Mr + k]) / hi[0]
                  - (preT[(i+1)*Mp*Mr + j*Mr + k] - preT[i*Mp*Mr + j*Mr + k]) / hi_1[0])
                  / (hi[0]/2 + hi_1[0]/2);

          Qxyz[1] = kxyz[1] * ((preT[i*Mp*Mr + j*Mr + k] - preT[i*Mp*Mr + (j-1)*Mr + k]) / hi[1]
                  - (preT[i*Mp*Mr + (j+1)*Mr + k] - preT[i*Mp*Mr + j*Mr + k]) / hi_1[1])
                  / (hi[1]/2 + hi_1[1]/2);

          Qxyz[2] = kxyz[2] * ((preT[i*Mp*Mr + j*Mr + k] - preT[i*Mp*Mr + j*Mr + (k-1)]) / hi[2]
                  - (preT[i*Mp*Mr + j*Mr + (k+1)] - preT[i*Mp*Mr + j*Mr + k]) / hi_1[2])
                  / (hi[2]/2 + hi_1[2]/2);
*/

          QQ[0] = - kxyz[0] * (preT[(i+1)*Mp*Mr + j*Mr + k] - preT[i*Mp*Mr + j*Mr + k]) / hi[0];
          QQ[1] = - kxyz[1] * (preT[i*Mp*Mr + (j+1)*Mr + k] - preT[i*Mp*Mr + j*Mr + k]) / hi[1];
          QQ[2] = - kxyz[2] * (preT[i*Mp*Mr + j*Mr + (k+1)] - preT[i*Mp*Mr + j*Mr + k]) / hi[2];

          Qxyz[0] = kxyz[0] * ( 2 / hi[0] / (hi_1[0]+hi[0]) * preT[(i+1)*Mp*Mr + j*Mr + k]
                  - 2 / hi_1[0] / hi[0] * preT[i*Mp*Mr + j*Mr + k]
                  + 2 / hi_1[0] / (hi_1[0] + hi[0]) * preT[(i-1)*Mp*Mr + j*Mr + k] );

          Qxyz[1] = kxyz[1] * ( 2 / hi[1] / (hi_1[1]+hi[1]) * preT[i*Mp*Mr + (j+1)*Mr + k]
                  - 2 / hi_1[1] / hi[1] * preT[i*Mp*Mr + j*Mr + k]
                  + 2 / hi_1[1] / (hi_1[1] + hi[1]) * preT[i*Mp*Mr + (j-1)*Mr + k] );

          Qxyz[2] = kxyz[2] * ( 2 / hi[2] / (hi_1[2]+hi[2]) * preT[i*Mp*Mr + j*Mr + (k+1)]
                  - 2 / hi_1[2] / hi[2] * preT[i*Mp*Mr + j*Mr + k]
                  + 2 / hi_1[2] / (hi_1[2] + hi[2]) * preT[i*Mp*Mr + j*Mr + (k-1)] );

          rho  = hc[i*Mp*Mr + j*Mr + k]->GetDensity();
          C    = hc[i*Mp*Mr + j*Mr + k]->GetCapacity();

          RadQ = hc[i*Mp*Mr + j*Mr + k]->GetHeatSource() * hc[i*Mp*Mr + j*Mr + k]->GetDensity();
          //RadQ = hc[i*Mp*Mr + j*Mr + k]->GetHeatSource();

          Q = (Qxyz[0] + Qxyz[1] + Qxyz[2] + RadQ) / rho / C;

          T[i*Mp*Mr + j*Mr + k] = preT[i*Mp*Mr + j*Mr + k] + dt * Q;

          hc[i*Mp*Mr + j*Mr + k]->SetHeat(Qxyz[0], Qxyz[1], Qxyz[2]);
          hc[i*Mp*Mr + j*Mr + k]->SetTime(time);
          hc[i*Mp*Mr + j*Mr + k]->SetPreTemperature(preT[i*Mp*Mr + j*Mr + k]);
          hc[i*Mp*Mr + j*Mr + k]->SetPostTemperature(T[i*Mp*Mr + j*Mr + k]);

          // *************************************************
          // DISPLAY:
          if (it%5000==1 && i==268 && j==2 && k==3)
            //tree->Fill();
            cout << fixed << setprecision(2) << time << "    "
                 << setw(5) << left << i
                 << setw(5) << left << j
                 << setw(5) << left << k
                 //<< Qxyz[0] << " " << Qxyz[1] << " " << Qxyz[2] << " "
                 << setw(10) << left << setprecision(4) << T[i*Mp*Mr + j*Mr + k] << endl;


          // *************************************************
          // RESULT:
          if (it%(20000*2)==1) {
            output << fixed << setprecision(2) << time << "    "; 
            output << setw(5) << left << i
                   << setw(4) << left << j
                   << setw(4) << left << k
                   << setw(12) << right << scientific << setprecision(4) << QQ[0]
                   << setw(12) << right << scientific << setprecision(4) << QQ[1]
                   << setw(12) << right << scientific << setprecision(4) << QQ[2]
                   << setw(12) << right << fixed << setprecision(4) << T[i*Mp*Mr + j*Mr + k] << "\n";
          }
        }
      }
    }

    // boundary condition
    for (int i=0; i<Mz; i++) {
      for (int j=0; j<Mp; j++) {
        T[i*Mp*Mr + j*Mr + 0 ] = T[i*Mp*Mr + j*Mr + 1];
        //T[i*Mp*Mr + j*Mr + (Mr-1)] = T[i*Mp*Mr + j*Mr + (Mr-2)];
        T[i*Mp*Mr + j*Mr + (Mr-1)] = 4.5;
      }
    }

    for (int i=0; i<Mz; i++) {
      for (int k=0; k<Mr; k++) {
        T[i*Mp*Mr + 0*Mr  + k] = T[i*Mp*Mr + (Mp-2)*Mr + k];
        T[i*Mp*Mr + (Mp-1)*Mr + k] = T[i*Mp*Mr + 1*Mr + k];
      }
    }

    for (int j=0; j<Mp; j++) {
      for (int k=0; k<Mr; k++) {
        T[(Mz-1)*Mp*Mr + j*Mr  + k] = T[(Mz-2)*Mp*Mr + j*Mr + k];
        T[0*Mp*Mr + j*Mr  + k] = T[1*Mp*Mr + j*Mr + k];
      }
    }

    for (int j=0; j<Mp; j++) {
      T[(Mz-1)*Mp*Mr + j*Mr + 1] = 4.5;
      T[0*Mp*Mr + j*Mr + 1] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 3] = 4.5;
      T[0*Mp*Mr + j*Mr + 3] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 5] = 4.5;
      T[0*Mp*Mr + j*Mr + 5] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 7] = 4.5;
      T[0*Mp*Mr + j*Mr + 7] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 9] = 4.5;
      T[0*Mp*Mr + j*Mr + 9] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 11] = 4.5;
      T[0*Mp*Mr + j*Mr + 11] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 13] = 4.5;
      T[0*Mp*Mr + j*Mr + 13] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 15] = 4.5;
      T[0*Mp*Mr + j*Mr + 15] = 4.5;
      T[(Mz-1)*Mp*Mr + j*Mr + 17] = 4.5;
      T[0*Mp*Mr + j*Mr + 17] = 4.5;
    }

  }

  output.close();
  //file->Write();
  //file->Close();
  
  return 0;
}
