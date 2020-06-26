#ifndef speckle_H
#define speckle_H

#include "std_include.h"
#include "noiseSource.h"

class speckle
    {
    public:
        
        void randomField(vector<vector<double> > &randField, double L, int nDisc, double sigma);
        double normal(double factor);
        void fft(int sign, vector<complex<double> > &zs);

    protected:

        bool Reproducible;
        noiseSource noise;

    };
#endif
