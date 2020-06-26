#include "std_include.h"

#include "Simulation.h"

#include "voronoiQuadraticEnergy.h"
#include "voronoiQuadraticEnergyNotConfluent.h"
#include "selfPropelledParticleDynamics.h"
#include "DatabaseNetCDFSPV.h"

#include <cstring>

using namespace std;


int main(){

  int last=1;
  int scriptsIni=1;
    for(int start=last; start<=last; start++){

      int v1, v2, v3, v4, v9, v12, v16, v17;
      double v5, v6, v7, v8, v10, v11, v13, v14, v15;
      //string access_dados="/home/diogo/MEGA/cenas/GPU/Results/scripts24/dados.txt";
      //string access_dados="/media/hdd/Results/scripts38/dados.txt";
      string access_dados="/media/hdd/Results_non_confluent/scripts"+to_string(scriptsIni)+"/dados.txt";
      //string access_dados="/media/hdd/Results/GPU/scripts1/dados.txt";
      //string access_dados="/home/diogo/MEGA/cenas/CellGPU_Substrate2/dados.txt";

      ifstream file(access_dados);

      for(int line=0; line<start; line++){
	file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11 >> v12 >> v13 >> v14 >> v15 >> v16 >> v17;
      }

      file.close();

      char data[1000000]; 
      //sprintf(data, "/home/diogo/MEGA/cenas/GPU/Results/scripts24/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/scripts38/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/GPU/scripts11/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/home/diogo/MEGA/cenas/CellGPU_Substrate2/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", v9, v10, v11, v6, v7, v8);
      sprintf(data, "/media/hdd/Results_non_confluent/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v17+1);

      SPVDatabaseNetCDF nc(v1,data,NcFile::ReadOnly);
      
        int N=v1;
        int a=N;
        int b=nc.GetNumRecs();
	
        int amostras=1;//v14;
        //int start=2;
        int fim=start;
	int time_total=b;
	int start_time=b-1;
        string filename;
	int val=1;
	cout<<"Size: "<<b<<endl;

	Dscalar boxsize=v14;
	Dscalar CellRadius=v15;
	EOMPtr spp = make_shared<selfPropelledParticleDynamics>(N);
        shared_ptr<VoronoiQuadraticEnergy> spv = make_shared<VoronoiQuadraticEnergyNotConfluent>(N, false, boxsize);
	spv->setCellRadiusUniform(CellRadius);
        spv->setCellAdaptationTimeUniform(v10);
	    
	//cout<<"here"<<endl;
	SimulationPtr sim = make_shared<Simulation>();
	sim->setConfiguration(spv);
	sim->addUpdater(spp,spv);
	bool initializeGPU=false;
	sim->setCPUOperation(!initializeGPU);
	
        FILE *f;
	    	        	
	for(int nn = start_time; nn<time_total; nn++){

	  cout<<"NEW TIME----------------"<<" "<<nn<<endl;
	    
	  for(int avg=(start-1)*amostras+1; avg<=(start-1)*amostras+1; avg++){
	    
            int curr_line=(avg-1)/amostras+1;           
            //int v1, v2, v3, v4, v9, v12, v14;
            //double v5, v6, v7, v8, v10, v11, v13;
            //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts20/dados.txt";
	    file.open(access_dados);

            for(int line=0; line<curr_line; line++){
	      file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11 >> v12 >> v13 >> v14 >> v15 >> v16 >> v17;
            }

            file.close();

            string add_JOB_num=to_string(avg);
            //if(avg<10)add_JOB_num="00"+to_string(avg);
            //else if(avg<100)add_JOB_num="0"+to_string(avg);
            //else add_JOB_num=to_string(avg);
           
            char dataname[1000000]; 
            //sprintf(dataname, "/home/diogo/MEGA/cenas/GPU/Results/scripts24/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/scripts38/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/GPU/scripts10/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
            //sprintf(dataname, "/home/diogo/MEGA/cenas/CellGPU_Substrate2/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", v9, v10, v11, v6, v7, v8);
	    sprintf(dataname, "/media/hdd/Results_non_confluent/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, avg);
	    
            SPVDatabaseNetCDF ncdat(N,dataname,NcFile::ReadOnly);
            //Check if the file exists in the output folder. if it does then do the scan
            if ((f = fopen(dataname, "r")) == NULL){
	      printf("Error! opening file\n");
	      printf("%s\n",dataname);
	      exit (911);
	      return -1;
            }
            else{
	      fclose(f);
            }

	    ncdat.ReadState(spv, nn, false);

	    filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/snap00"+to_string(val)+".pov";
	    //filename="/home/diogo/MEGA/cenas/CellGPU_Substrate2/snap00"+to_string(val)+".pov";
	    val++;
	    //filename="/media/hdd/PovRay/2048/snap00"+to_string(val)+".pov";
	    ofstream file1;
	    file1.open (filename);
	    spv->writeTriangulation(file1);
	    file1.close();
	  }
	}
    }
    return 0;
}
