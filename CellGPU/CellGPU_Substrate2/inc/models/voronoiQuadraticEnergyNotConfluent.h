#ifndef VoronoiQuadraticEnergyNotConfluent_H
#define VoronoiQuadraticEnergyNotConfluent_H

#include "voronoiQuadraticEnergy.h"

/*! \file voronoiQuadraticEnergy.h */
//!Implement a 2D Voronoi model, with and without some extra bells and whistles, using kernels in \ref spvKernels
/*!
 *A child class of voronoiModelBase, this implements a Voronoi model in 2D. This involves mostly calculating
  the forces in the Voronoi model. Optimizing these procedures for
  hybrid CPU/GPU-based computation involves declaring and maintaining several related auxiliary
  data structures that capture different features of the local topology and local geoemetry for each
  cell.
 */
class VoronoiQuadraticEnergyNotConfluent : public VoronoiQuadraticEnergy
    {
    public:
        //!initialize with random positions in a square box
        VoronoiQuadraticEnergyNotConfluent(int n,bool reprod, Dscalar boxsize);

        //!Blank constructor
        VoronoiQuadraticEnergyNotConfluent(){};

        //!Initialize voronoiQuadraticEnergy and call the initializer chain
        void initializeVoronoiQuadraticEnergyNotConfluent(int n);

        //!compute the geometry and get the forces
        virtual void computeForces();

        //!compute the quadratic energy functional
        virtual Dscalar computeEnergy();

        //cell-dynamics related functions...these call functions in the next section
        //in general, these functions are the common calls, and test flags to know whether to call specific versions of specialty functions

        //CPU functions
        //!Compute the net force on particle i on the CPU
        virtual void computeVoronoiNonConfluentForceCPU(int i);
	virtual void computeVoronoiNonConfluentForceCPU(int i, int j, Dscalar2 &fij);

	virtual void AreaPeriCalc_for_Forces(Dscalar2 &dPm, Dscalar2 &dAm, int type1, int type2, Dscalar2 vcur, Dscalar2 vnext, Dscalar l, bool sign, Dscalar teta);

	void computeLocalStressTensor(int i, Matrix2x2 &sint);
	void computeSwimStressTensor(int i, Matrix2x2 &sswim);

	virtual void reportError(ofstream &outfile);
    protected:

    //be friends with the associated Database class so it can access data to store or read
    friend class SPVDatabaseNetCDF;
    };

#endif
