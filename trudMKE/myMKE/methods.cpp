//
// Created by Artem Burcha on 12.11.2017.
//

#include "classdefinition.h"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

GlobalMatrix::GlobalMatrix() {cout<<"NO PARAMETERS INPUT IN CONSTRUCTOR!!!"; exit(0);}

GlobalMatrix::GlobalMatrix(double aa, double bb, double cc, double dd,
                           double inalpha, double inbeta, int elem, double L, double a):
a(aa),b(bb),c(cc),d(dd),alpha(inalpha),beta(inbeta),nodeAmount(elem+1),elementLength(L),special(a)
{
    GM = new double* [nodeAmount];
    for (int i = 0; i < nodeAmount; i++) {
        GM[i] = new double[nodeAmount + 1];
    }
    InitMatrix();
}

void GlobalMatrix::InitMatrix() {
        for (int i = 0; i < nodeAmount; i++) {
            for (int j = 0; j < nodeAmount +1 ; j++)
                GM[i][j] = 0;
        }
}

void GlobalMatrix::PrintMatrix() {
    cout.precision(3);
    cout.setf(ios::left);
    cout.width(7);
    for (int i = 0; i < nodeAmount; i++) {
        for (int j = 0; j < nodeAmount +1 ; j++)
            cout<<GM[i][j]<<"\t";
        cout<<endl;
    }
}

void GlobalMatrix::Ensembling() {
   // double x=2.0;
    for (int i=0; i<nodeAmount-1; i++){

        GM[i][i]     += a;
        GM[i][i+1]   += b;
        GM[i+1][i]   += c;
        GM[i+1][i+1] += d;
        GM[i][nodeAmount]   += alpha;
        GM[i+1][nodeAmount] += beta;
    }

    GM[0][0] = -3.0;
    GM[1][0] = 0.0;
    GM[nodeAmount-2][nodeAmount-1]=0.0;
    GM[nodeAmount-1][nodeAmount-1]=3.0;

    //GM[0][nodeAmount]-=a*0; //u(2)=0; GM[1][nodeAmount] doesn't change too
    GM[nodeAmount-2][nodeAmount] -=b*10.0; //u(15) = 10;
    GM[nodeAmount-1][nodeAmount] -=d*10.0;
}

void GlobalMatrix::Gauss() {
        int      i, j, k, maxInd;
        double    max;
        double    tmp;
        double    eps = 0.00001;
        /*Direct*/
        for(i = 0; i < nodeAmount; i++) {
            max = fabs(GM[i][i]);
            maxInd = i;
            // Search max diag element
            for(k = i + 1; k < nodeAmount; k++) {
                if(fabs(GM[k][i]) > max) {
                    maxInd = k;
                    max = fabs(GM[k][i]);
                }
            }
            // if max is zero
            if(max < eps) {
                puts("Can'GM solve system of equations: zero column");
                exit(5);
            }
            // Swap current and max lines
            for(j = 0; j < nodeAmount + 1; j++) {
                tmp = GM[maxInd][j];
                GM[maxInd][j] = GM[i][j];
                GM[i][j] = tmp;
            }
            // Divide each element in diag line
            tmp = GM[i][i];
            for(j = i; j < nodeAmount + 1; j++) {
                GM[i][j] /= tmp;
            }
            // Subtraction multiple first line from each line
            for(k = i + 1; k < nodeAmount; k++) {
                tmp = GM[k][i];
                for(j = i; j < nodeAmount + 1; j++) {
                    GM[k][j] -= GM[i][j] * tmp;
                }
            }
        }
        // Reverse
        for(i = nodeAmount-1; i >=0; i--) {
            double tmp = GM[i][nodeAmount];
            solution.push_back(tmp);
            for(j = 0; j < i; j++) {
                GM[j][nodeAmount] -= GM[j][i] * tmp;
            }
        }
}

double GlobalMatrix::analyticSolution(double x) {
    double u;
    u = (-1.0/(exp(26.0*sqrt(15.0)/3.0) - 1.0)) *
            (2.0*(-4.0*exp((x-2.0)*sqrt(15.0))
                - exp(-19.0*sqrt(15.0)/3.0 + sqrt(15.0)*x)
                - exp(3.0*sqrt(15.0)+2.0/3.0*sqrt(15.0)*x)
                + exp((2.0*x - 17.0)*sqrt(15.0)/3.0)
                + exp((11.0+x)*sqrt(15.0)/3.0)
                + 4.0*exp((x-2.0)*sqrt(15.0)/3.0))
             *exp((-2.0*x+17.0)*sqrt(15.0)/3.0));
    return u;
}

void GlobalMatrix::Compare() {
    cout.precision(6);
    //cout.setf(ios::fixed);
    ofstream graph("../data.txt");
    graph<<2.0<<"\t"<<0.0<<endl;
    double currentX = 2.0;
    double maxDiff = 0.0;
    double difference;
    double currentAnalytic;
    cout<<"result"<<"\t\t"<<"analytic"<<"\t\t"<<"difference"<<endl;
    for (int i=2; i<nodeAmount; i++)
    {
        currentX += elementLength;
        currentAnalytic = analyticSolution(currentX);
        difference = fabs(currentAnalytic - solution[nodeAmount-i]);
        cout<<solution[nodeAmount-i]<<"\t\t"<<currentAnalytic<<"\t\t"<<difference<<endl;
        maxDiff = (maxDiff < difference) ? difference : maxDiff;
        graph<<currentX<<"\t"<<solution[nodeAmount-i]<<endl;
    }
    graph<<15.0<<"\t"<<10.0<<endl;
    cout<<endl<<"Max difference = "<<maxDiff<<endl;
}
