#include "std_include.h"

#include "cuda_runtime.h"
#include "cuda_profiler_api.h"


#include "Simulation.h"
#include "voronoiQuadraticEnergy.h"
#include "voronoiQuadraticEnergyNotConfluent.h"
#include "selfPropelledParticleDynamics.h"
#include "DatabaseNetCDFSPV.h"
#include "speckle.h"

/*!
This file compiles to produce an executable that can be used to reproduce the timing information
in the main cellGPU paper. It sets up a simulation that takes control of a voronoi model and a simple
model of active motility
*/
int main(int argc, char*argv[])
{
    //...some default parameters
    int numpts = 120; //number of cells
    int USE_GPU = -1; //0 or greater uses a gpu, any negative number runs on the cpu
    int c;
    int tSteps = 10000; //number of time steps to run after initialization
    int initSteps = 100; //number of initialization steps

    Dscalar dt = 0.001; //the time step size
    Dscalar p0 = 2.5;  //the preferred perimeter
    Dscalar a0 = 4.0;  // the preferred area
    Dscalar v0 = 0.1;  // the self-propulsion
    Dscalar Dr = 1.0;  // the rotational diffusion
    int dx=1;
    int program_switch = -1; //various settings control output
    int SubIntType = 0;
    Dscalar IncrementationValue=0.0; 
    bool rnd_Sub=false;
    Dscalar tau=1;
    Dscalar boxsize=10;
    Dscalar CellRadius=1;//if this value is smaller or equal to 0 then use the standard VoronoiQuadraticEnergy branch, otherwise use the NonConfluent branch
    int CellClusterSize=numpts;

    double sigmac = 0.;
    char fn[256];
    sprintf(fn,"snap001.pov");
    ofstream output(fn);
    
    //The defaults can be overridden from the command line
    while((c=getopt(argc,argv,"a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:")) != -1)
        switch(c)
        {
            case 'a': numpts = atoi(optarg); break;
            case 'b': tSteps = atoi(optarg); break;
            case 'c': USE_GPU = atoi(optarg); break;
            case 'd': initSteps = atoi(optarg); break;
            case 'e': dt = atof(optarg); break;
            case 'f': p0 = atof(optarg); break;
            case 'g': a0 = atof(optarg); break;
            case 'h': v0 = atof(optarg); break;
            case 'i': SubIntType = atoi(optarg); break;
            case 'j': tau = atof(optarg); break;
            case 'k': IncrementationValue = atof(optarg); break;
            case 'l': dx = atoi(optarg); break;
            case 'm': sigmac = atof(optarg); break;
            case 'n': boxsize = atof(optarg); break;
            case 'o': CellRadius = atof(optarg); break;
            case 'p': CellClusterSize = atoi(optarg); break;
            case '?':
                    if(optopt=='c')
                        std::cerr<<"Option -" << optopt << "requires an argument.\n";
                    else if(isprint(optopt))
                        std::cerr<<"Unknown option '-" << optopt << "'.\n";
                    else
                        std::cerr << "Unknown option character.\n";
                    return 1;
            default:
                       abort();
        };

    clock_t t1,t2; //clocks for timing information
    bool reproducible = true; // if you want random numbers with a more random seed each run, set this to false
    //check to see if we should run on a GPU
    bool initializeGPU = true;
    if (USE_GPU >= 0)
        {
	  initializeGPU = false;
        //bool gpu = chooseGPU(USE_GPU);
        //if (!gpu) return 0;
        //cudaSetDevice(USE_GPU);
        }
    else
        initializeGPU = false;

    //possibly save output in netCDF format
    char dataname[256];
    sprintf(dataname,"./test_voronoi.nc");
    int Nvert = numpts;
    SPVDatabaseNetCDF ncdat(Nvert,dataname,NcFile::Replace);

    //define an equation of motion object...here for self-propelled cells
    EOMPtr spp = make_shared<selfPropelledParticleDynamics>(numpts);
    
    //define a voronoi configuration with a quadratic energy functional
    //shared_ptr<VoronoiQuadraticEnergyNotConfluent> spv  = make_shared<VoronoiQuadraticEnergyNotConfluent>(numpts, reproducible, boxsize);
    //shared_ptr<VoronoiQuadraticEnergy> spv  = make_shared<VoronoiQuadraticEnergy>(numpts, reproducible);
    
    shared_ptr<VoronoiQuadraticEnergy> spv;
    if(CellRadius>0) spv  = make_shared<VoronoiQuadraticEnergyNotConfluent>(numpts, reproducible, boxsize);
    else spv  = make_shared<VoronoiQuadraticEnergy>(numpts, reproducible);
    
    //set the cell preferences to uniformly have A_0 = 1, P_0 = p_0
    spv->setCellPreferencesUniform(a0,p0);
    //set the cell activity to have D_r = 1. and a given v_0
    spv->setv0Dr(v0,1.0);
    //set the substrate preferences values
    spv->setSubstratePreferences(0, 0.0, 1, 0.0);
    //set the cell moduli preferences to uniformly have a given K_a and K_p
    spv->setModuliUniform(1.0, 1.0);
    //set the cell radius to uniform values
    spv->setCellRadiusUniform(CellRadius);
    //set the cell adaptation time to uniformly have tau = tau
    spv->setCellAdaptationTimeUniform(tau);
    //enable simulation of small cell clusters
    if(CellClusterSize<numpts)spv->setExclusionsRadius(CellClusterSize);

    //combine the equation of motion and the cell configuration in a "Simulation"
    SimulationPtr sim = make_shared<Simulation>();
    sim->setConfiguration(spv);
    sim->addUpdater(spp,spv);
    //set the time step size
    sim->setIntegrationTimestep(dt);
    //initialize Hilbert-curve sorting... can be turned off by commenting out this line or seting the argument to a negative number
    //sim->setSortPeriod(initSteps/10);
    //set appropriate CPU and GPU flags
    sim->setCPUOperation(!initializeGPU);
    sim->setReproducible(reproducible);

    ncdat.WriteState(spv);
    //spv->writeTriangulation(output);
    //exit (911);
    //run for a few initialization timesteps
    //printf("starting initialization\n");
    for(int ii = 0; ii < initSteps; ++ii)
        {
        sim->performTimestep();
	ncdat.WriteState(spv);
	//spv->writeTriangulation(output);
	//exit (911);
        };
    printf("Finished with initialization\n");
    cout << "current q = " << spv->reportq() << endl;
    //the reporting of the force should yield a number that is numerically close to zero.
    spv->reportMeanCellForce(true);
    double arguments[dx*dx];
    //spv->writeTriangulation(output);
    //exit (911);

    if(rnd_Sub==false)spv->setSubstratePreferences(SubIntType, IncrementationValue, dx, sigmac);
    else{ 
        vector<vector<double> > randField(dx, vector<double>(dx,0));
        speckle rndfield;
        rndfield.randomField(randField, sqrt(numpts), dx, sigmac);

        for(int u=0; u<dx; u++){
            for(int k=0; k<dx; k++){
                int calc_site=k+u*dx;
                arguments[calc_site]=randField[k][u];
            }
        }
        spv->setSubstratePreferencesRnd(SubIntType, IncrementationValue, dx, arguments, sigmac);
    }

    switch(SubIntType)
    {
    case 1:
        spv->setCellPreferencesSubstrate(a0, p0, true, dx);
    break;
    case 2:
        spv->setCellPreferencesSubstrate(a0, p0, false, dx);
    break;
    case 3:
        spv->setv0DrSubstrate(v0, 1.0, dx);
    break;
    }
  
    //run for additional timesteps, and record timing information
    t1=clock();
    int mult=10;
    int val=1;
    int value_calc=10;
    int count_print=1000;
    int count_print2=100000;
    for(int ii = 0; ii < tSteps; ++ii)
        {

        if(ii%1000 ==0)
            {
            printf("timestep %i\t\t energy %f \n",ii,spv->computeEnergy());
            };
        sim->performTimestep();
	ncdat.WriteState(spv);
        if(program_switch <0)
	  {
	      
	    int power=pow(10, val);
            if(ii%power==0){
              switch(val){
                case 2:
                  mult=100;
                break;
                case 3:
                  mult=200;
                break;
                case 6:
                  mult=300;
                break;
              }

            value_calc=(int)(pow(10, val+1)/mult);
            val++;
            }

            if(ii%value_calc==0)ncdat.WriteState(spv);     
            };
        };
    cout<<"energy: "<<spv->computeEnergy()<<endl;
    t2=clock();
    Dscalar steptime = (t2-t1)/(Dscalar)CLOCKS_PER_SEC/tSteps;
    cout << "timestep ~ " << steptime << " per frame; " << endl;
    cout << spv->reportq() << endl;
    spv->reportMeanCellForce(true);
    spv->writeTriangulation(output);

/*
    if(initializeGPU)
        cudaDeviceReset();*/
    return 0;
};
