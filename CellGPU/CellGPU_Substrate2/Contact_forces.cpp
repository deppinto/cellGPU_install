#include "std_include.h"

#include "Simulation.h"

#include "voronoiQuadraticEnergy.h"
#include "voronoiQuadraticEnergyNotConfluent.h"
#include "selfPropelledParticleDynamics.h"
#include "DatabaseNetCDFSPV.h"
#include "Matrix.h"

#include <cstring>

using namespace std;

#define PI  3.141592653589793

int main(){

  int last=1;
  int scriptsIni=4;
  int level=0;
    for(int start=1; start<=last; start++){

      int v1, v2, v3, v4, v9, v12, v16, v17;
      double v5, v6, v7, v8, v10, v11, v13, v14, v15;
      //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts31/dados.txt";
      //string access_dados="/media/hdd/Results/scripts42/dados.txt";
      string access_dados="/media/hdd/Results_non_confluent/scripts"+to_string(scriptsIni)+"/dados.txt";
      //string access_dados="/media/hdd/Results/GPU/scripts23/dados.txt";
      //string access_dados="/home/diogo/MEGA/cenas/GPU/Results/scripts24/dados.txt";

      ifstream file(access_dados);

      for(int line=0; line<start; line++){
	file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11 >> v12 >> v13 >> v14 >> v15 >> v16 >> v17;
      }

      file.close();

      char data[1000000]; 
      //sprintf(dataname, "/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/scripts42/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      sprintf(data, "/media/hdd/Results_non_confluent/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v17+1);
      //sprintf(data, "/media/hdd/Results/GPU/scripts23/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/home/diogo/MEGA/cenas/GPU/Results/scripts24/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
 
      SPVDatabaseNetCDF nc(v1,data,NcFile::ReadOnly);
      
        int N=v1;
        int a=N;
        int b=nc.GetNumRecs();

	double KA=1.0;
	double KP=1.0;
        int amostras=v17;
        int fim=start;
	int time_total=b;
	int start_time=0;
        vector <double> time(b,0);
	vector <double> AvgNormalStress(b,0);
	vector <double> AvgShearStress(b,0);
	vector <double> AvgForcex(b,0);
	vector <double> AvgForcey(b,0);
	vector <double> AvgForce(b,0);
	vector <double> AvgEnergy(b,0);
        string filename;
	Dscalar boxsize=v14;
	Dscalar CellRadius=v15;

	EOMPtr spp = make_shared<selfPropelledParticleDynamics>(N);
        shared_ptr<VoronoiQuadraticEnergyNotConfluent> spv  = make_shared<VoronoiQuadraticEnergyNotConfluent>(N, false, boxsize);
        spv->setCellRadiusUniform(CellRadius);
        spv->setCellAdaptationTimeUniform(v10);
	SimulationPtr sim = make_shared<Simulation>();
	sim->setConfiguration(spv);
	sim->addUpdater(spp,spv);
	bool initializeGPU=false;
	sim->setCPUOperation(!initializeGPU);

        FILE *f;

	Dscalar minSS=100000000;
        Dscalar maxSS=-100000000;
        Dscalar minNS=100000000;
        Dscalar maxNS=-100000000;
        Dscalar minFX=100000000;
        Dscalar maxFX=-100000000;
        Dscalar minFY=100000000;
        Dscalar maxFY=-100000000;
        Dscalar minF=100000000;
        Dscalar maxF=-100000000;
        Dscalar minE=100000000;
        Dscalar maxE=-100000000;
	bool Second=false;
        vector <double> HistNormalStress;
        vector <double> HistShearStress;
        vector <double> HistForcex;
        vector <double> HistForcey;
	vector <double> HistForce;
        vector <double> HistEnergy;
        Dscalar binSS=0;
        Dscalar binNS=0;
        Dscalar binFX=0;
        Dscalar binFY=0;
	Dscalar binF=0;
        Dscalar binE=0;
	int numBins=100;
	ofstream file1;
	
	for(int nn = start_time; nn<time_total; nn++){

	  cout<<"NEW TIME----------------"<<" "<<nn <<endl;
	    
	  for(int avg=(start-1)*amostras+1; avg<=fim*amostras; avg++){

		  //cout<<"NEW TIME----------------"<<" "<<avg<<endl;
		  if(Second==true)
		  {
			  binE=abs(maxE-minE)/numBins;
			  HistEnergy.resize(numBins, 0);

                          binFX=abs(maxFX-minFX)/numBins;
                          HistForcex.resize(numBins, 0);

			  binFY=abs(maxFY-minFY)/numBins;
                          HistForcey.resize(numBins, 0);

			  binF=abs(maxF-minF)/numBins;
                          HistForce.resize(numBins, 0);

			  binSS=abs(maxSS-minSS)/numBins;
                          HistShearStress.resize(numBins, 0);

                          binNS=abs(maxNS-minNS)/numBins;
                          HistNormalStress.resize(numBins, 0);
		  }

            int curr_line=(avg-1)/amostras+1;           
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
            //sprintf(dataname, "/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/scripts42/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    sprintf(dataname, "/media/hdd/Results_non_confluent/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, avg);
	    //sprintf(dataname, "/media/hdd/Results/GPU/scripts23/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
            //sprintf(dataname, "/home/diogo/MEGA/cenas/GPU/Results/scripts24/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
 
	    
            SPVDatabaseNetCDF ncdat(N,dataname,NcFile::ReadOnly);
            //Check if the file exists in the output folder. if it does then do the scan
            if ((f = fopen(dataname, "r")) == NULL){
	      printf("Error! opening file: %d\n", avg);
	      printf("%s\n",dataname);
	      throw exception();
	      return -1;
            }
            else{
	      fclose(f);
            }

	    ncdat.ReadState(spv, nn, true);

	    //spv->globalTriangulationCGAL();
	    //spv->ComputeGeometryNonConfluent();

	    //read in all the data we'll need
            ArrayHandle<int> h_vn(spv->voroNum,access_location::host,access_mode::read);
            ArrayHandle<Dscalar2> h_AP(spv->AreaPeri,access_location::host,access_mode::read);
	    ArrayHandle<Dscalar2> h_APP(spv->AreaPeriPreferences,access_location::host,access_mode::read);

	    Dscalar2 fij;
	    Dscalar2 forceTotal;
	    Dscalar2 meanForce;
	    meanForce.x=0;
	    meanForce.y=0;
	    Dscalar Energy;
	    Dscalar meanEnergy=0;
	    forceTotal.x=0;
	    forceTotal.y=0;
	    Dscalar FTotal=0;
	    Dscalar meanF=0;
	    Matrix2x2 LocalStressTensor;
	    Matrix2x2 SwimStressTensor;

	    Matrix2x2 StressTensor;
	    Dscalar normalStress;
	    Dscalar shearStress;

            Matrix2x2 meanStressTensor;
	    meanStressTensor=0;
            Dscalar meanNormalStress=0;
            Dscalar meanShearStress=0;

	    int binVal;
	    Dscalar totalArea=0;
	    vector<double> colorE(N, 0);
            vector<double> colorFX(N, 0);
            vector<double> colorFY(N, 0);
            vector<double> colorF(N, 0);
            vector<double> colorSS(N, 0);
            vector<double> colorNS(N, 0);

        for(int i=0; i< N; i++)
	{
		forceTotal.x=0;
		forceTotal.y=0;
		Energy=0;
		FTotal=0;
		int neighMax=h_vn.data[i];
		for(int j=0; j<neighMax; j++)
		{
			fij.x=0;
			fij.y=0;
			spv->computeVoronoiNonConfluentForceCPU(i, j, fij);
			forceTotal.x+=fij.x;
			forceTotal.y+=fij.y;
		}
		FTotal=sqrt(forceTotal.x*forceTotal.x+forceTotal.y*forceTotal.y);
		meanF+=FTotal;
		meanForce.x+=forceTotal.x;
		meanForce.y+=forceTotal.y;

		Energy+= KA * (h_AP.data[i].x-h_APP.data[i].x)*(h_AP.data[i].x-h_APP.data[i].x);
		Energy+= KP * (h_AP.data[i].y-h_APP.data[i].y)*(h_AP.data[i].y-h_APP.data[i].y);
	        meanEnergy += Energy;

		spv->computeLocalStressTensor(i, LocalStressTensor);
		spv->computeSwimStressTensor(i, SwimStressTensor);

		StressTensor=LocalStressTensor+SwimStressTensor;
		normalStress=0.5*(StressTensor.x11+StressTensor.x22);
		shearStress=0.5*(StressTensor.x21+StressTensor.x12);

		Dscalar cellArea=h_AP.data[i].x;
                meanStressTensor+=cellArea*LocalStressTensor+cellArea*SwimStressTensor;
		totalArea+=cellArea;

		if(nn==time_total-1 && Second==false)
		{
        		if(Energy>maxE)maxE=Energy;
		        if(Energy<minE)minE=Energy;

		        if(forceTotal.x>maxFX)maxFX=forceTotal.x;
		        if(forceTotal.x<minFX)minFX=forceTotal.x;

		        if(forceTotal.y>maxFY)maxFY=forceTotal.y;
		        if(forceTotal.y<minFY)minFY=forceTotal.y;

                        if(FTotal>maxF)maxF=FTotal;
                        if(FTotal<minF)minF=FTotal;

	        	if(normalStress>maxNS)maxNS=normalStress;
	        	if(normalStress<minNS)minNS=normalStress;

	        	if(shearStress>maxSS)maxSS=shearStress;
	        	if(shearStress<minSS)minSS=shearStress;

			if(avg==fim*amostras && i==N-1){avg=(start-1)*amostras;Second=true;}
		}
		else if(nn==time_total-1 && Second==true)
		{
			  colorE[i]=Energy;
                	  colorFX[i]=forceTotal.x;
                	  colorFY[i]=forceTotal.y;
                	  colorF[i]=FTotal;
                	  colorSS[i]=shearStress;
                	  colorNS[i]=normalStress;

	 		  if(Energy>maxE)colorE[i]=maxE;
			  if(Energy<minE)colorE[i]=minE;
                          binVal=(colorE[i]-minE)/binE;
                          HistEnergy[binVal]++;

			  if(forceTotal.x>maxFX)colorFX[i]=maxFX;
                          if(forceTotal.x<minFX)colorFX[i]=minFX;
                          binVal=(colorFX[i]-minFX)/binFX;
                          HistForcex[binVal]++;

			  if(forceTotal.y>maxFY)colorFY[i]=maxFY;
                          if(forceTotal.y<minFY)colorFY[i]=minFY;
			  binVal=(colorFY[i]-minFY)/binFY;
                          HistForcey[binVal]++;

			  if(FTotal>maxF)colorF[i]=maxF;
                          if(FTotal<minF)colorF[i]=minF;
			  binVal=(colorF[i]-minF)/binF;
                          HistForce[binVal]++;

			  if(shearStress>maxSS)colorSS[i]=maxSS;
                          if(shearStress<minSS)colorSS[i]=minSS;
		  	  binVal=(colorSS[i]-minSS)/binSS;
                       	  HistShearStress[binVal]++;

			  if(normalStress>maxNS)colorNS[i]=maxNS;
                          if(normalStress<minNS)colorNS[i]=minNS;
		  	  binVal=(colorNS[i]-minNS)/binNS;
                       	  HistNormalStress[binVal]++;

			  if(avg==(start-1)*amostras+1 && i==N-1)
			  {
		          	filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/Esnap00"+to_string(start)+".pov";
			        file1.open (filename);
        			spv->writeTriangulation(file1, binE*numBins, minE, colorE, false);
        			file1.close();

                                filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/FXsnap00"+to_string(start)+".pov";
                                file1.open (filename);
                                spv->writeTriangulation(file1, binFX*numBins, minFX, colorFX, false);
                                file1.close();

                                filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/FYsnap00"+to_string(start)+".pov";
                                file1.open (filename);
                                spv->writeTriangulation(file1, binFY*numBins, minFY, colorFY, false);
                                file1.close();

                                filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/Fsnap00"+to_string(start)+".pov";
                                file1.open (filename);
                                spv->writeTriangulation(file1, binF*numBins, minF, colorF, false);
                                file1.close();

                                filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/SSsnap00"+to_string(start)+".pov";
                                file1.open (filename);
                                spv->writeTriangulation(file1, binSS*numBins, minSS, colorSS, false);
                                file1.close();

                                filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/NSsnap00"+to_string(start)+".pov";
                                file1.open (filename);
                                spv->writeTriangulation(file1, binNS*numBins, minNS, colorNS, false);
                                file1.close();
			  }
		}
	}
	meanStressTensor*=(1/totalArea);
        meanNormalStress=0.5*(meanStressTensor.x11+meanStressTensor.x22);
        meanShearStress=0.5*(meanStressTensor.x21+meanStressTensor.x12);

	AvgNormalStress[nn]+=meanNormalStress;
	AvgShearStress[nn]+=meanShearStress;

	AvgForcex[nn]+=meanForce.x/N;
	AvgForcey[nn]+=meanForce.y/N;
	AvgForce[nn]+=meanF/N;

	AvgEnergy[nn]=meanEnergy/N;

	//cout<<"STRESS: "<<meanNormalStress<<" "<<meanShearStress<<" "<<meanForce.x/N<<" "<<meanForce.y/N<<" "<<meanEnergy/N<<endl;

	time[nn]= spv->currentTime;
	}
	}	  


	filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/CF_"+to_string(start)+"_"+to_string(level)+".txt";
	//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Stripe/s35/Z_"+to_string(start)+".txt";
	file1.open (filename);
	for(int t_step=start_time; t_step<b; t_step++){
		
	file1 << time[t_step] << " " << AvgEnergy[t_step]/(double)amostras<<" "<< AvgForcex[t_step]/(double)amostras<<" "<< AvgForcey[t_step]/(double)amostras  <<" "<< AvgForce[t_step]/(double)amostras <<" "<< AvgShearStress[t_step]/(double)amostras <<" "<< AvgNormalStress[t_step]/(double)amostras <<endl;
	}
	file1.close();

        filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/Dists_"+to_string(start)+"_"+to_string(level)+".txt";
        file1.open (filename);
        for(int t_step=0; t_step<numBins; t_step++){

        file1 << minE+t_step*binE << " " << HistEnergy[t_step]/(amostras*N)<<" "<< minFX+t_step*binFX <<" "<< HistForcex[t_step]/(amostras*N)  <<" "<< minFY+t_step*binFY <<" "<< HistForcey[t_step]/(amostras*N)  <<" "<< minF+t_step*binF <<" "<< HistForce[t_step]/(amostras*N)  <<" "<< minSS+t_step*binSS <<" "<< HistShearStress[t_step]/(amostras*N)  <<" "<< minNS+t_step*binNS <<" "<< HistNormalStress[t_step]/(amostras*N)  <<endl;
        }
        file1.close();
    }
    return 0;
}
