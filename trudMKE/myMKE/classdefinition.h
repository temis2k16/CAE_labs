//
// Created by Artem Burcha on 11.11.2017.
//

#include <vector>

using namespace std;

class GlobalMatrix{
private:
    double a;
    double b;
    double c;
    double d;
    double alpha;
    double beta;
    int nodeAmount;
    double elementLength;
    double **GM;
    vector<double> solution;
    double special;

public:
    GlobalMatrix();
    GlobalMatrix(double aa, double bb, double cc, double dd, double inalpha, double inbeta,int elementAmount,double elementLength,double a);
    ~GlobalMatrix();
    void InitMatrix();
    void Ensembling();
    void Gauss();
    double analyticSolution(double x);
    void Compare();
    void PrintMatrix();
    void toFile();
};