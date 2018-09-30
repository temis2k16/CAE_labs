#include <iostream>
#include <cmath>
#include "classdefinition.h"
#define LEFTBORDER 2.0
#define RIGHTBORDER 15.0
#define N 40.0
//#define QUBIC  //LINEAR, QUBIC
#define LINEAR

using namespace std;

int main() {

    double A;
    double B;
    double C;
    double D;
    double ALPH;
    double BET;

#ifdef LINEAR
    //LL = Linear length of element
    double start = 6.0;
    double end = 10.0;
    double special = end - start;
    double LL = (RIGHTBORDER - LEFTBORDER)/(N);
    A = -5.0*LL/3.0 - 3.0/LL;
    B = 3.0/LL - 5.0*LL/6.0;
    C = 3.0/LL - 5.0*LL/6.0;
    D = -5.0*LL/3.0 - 3.0/LL;
    ALPH = -5.0*LL;
    BET = -5.0*LL;
    GlobalMatrix *linear = new GlobalMatrix(A,B,C,D,ALPH,BET,N,LL,special);
    linear->Ensembling();
    linear->PrintMatrix();
    linear->Gauss();
    cout<<endl;
    linear->PrintMatrix();
    cout<<LL<<endl;
    cout<<endl;
    linear->Compare();
#endif

#ifdef QUBIC
    double LQ = (RIGHTBORDER - LEFTBORDER)/N;
     A = - (5.0/3.0 * pow(LQ,6.0) + 135.0*pow(LQ,4.0) + 1728.0*pow(LQ,2.0) + 2268.0) /
             (LQ*(5.0*pow(LQ,2.0)+126.0)*(pow(LQ,2.0)+6.0));
     B = - (5.0*pow(LQ,6.0) - 90.0*pow(LQ,4.0) + 1944.0*pow(LQ,2.0) - 27216.0) /
             (12.0*LQ*(5.0*pow(LQ,2.0)+126.0)*(pow(LQ,2.0)+6.0));
     C =  - (5.0*pow(LQ,6.0) - 90.0*pow(LQ,4.0) + 1944.0*pow(LQ,2.0) - 27216.0) /
          (12.0*LQ*(5.0*pow(LQ,2.0)+126.0)*(pow(LQ,2.0)+6.0));
     D = - (5.0/3.0 * pow(LQ,6.0) + 135.0*pow(LQ,4.0) + 1728.0*pow(LQ,2.0) + 2268.0) /
         (LQ*(5.0*pow(LQ,2.0)+126.0)*(pow(LQ,2.0)+6.0));
     ALPH = -5.0*LQ*(pow(LQ,2.0) + 36.0)/(6.0*(pow(LQ,2.0)+6.0));
     BET = -5.0*LQ*(pow(LQ,2.0) + 36.0)/(6.0*(pow(LQ,2.0)+6.0));

    GlobalMatrix *quadratic = new GlobalMatrix(A,B,C,D,ALPH,BET,N,LQ,1);
    quadratic->Ensembling();
    cout<<LQ<<endl;
    quadratic->Gauss();
    cout<<endl;
    quadratic->Compare();
#endif
    return 0;
}