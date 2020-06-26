#include "std_include.h"

#include "Simulation.h"

#include "voronoiQuadraticEnergy.h"
#include "voronoiQuadraticEnergyNotConfluent.h"
#include "selfPropelledParticleDynamics.h"
#include "DatabaseNetCDFSPV.h"
#include "Matrix.h"

#include <cstring>

using namespace std;

int main(){

    int cont=0;
    int last=1;
    int scriptsIni=4;
    for(int start=1; start<=last; start++){

      int v1, v2, v3, v4, v9, v12, v16, v17;
      double v5, v6, v7, v8, v10, v11, v13, v14, v15;
      //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts35/dados.txt";
      //string access_dados="/media/hdd/Results/scripts39/dados.txt";
      //string access_dados="/media/hdd/Results/GPU/scripts5/dados.txt";
      string access_dados="/media/hdd/Results_non_confluent/scripts"+to_string(scriptsIni)+"/dados.txt";
      ifstream file(access_dados);

      for(int line=0; line<start; line++){
	file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11 >> v12 >> v13 >> v14 >> v15 >> v16 >> v17;
      }

      file.close();

      char data[1000000]; 
      //sprintf(dataname, "/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/scripts39/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/GPU/scripts5/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      sprintf(data, "/media/hdd/Results_non_confluent/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v17+1);
 
      SPVDatabaseNetCDF nc(v1,data,NcFile::ReadOnly);

        //int ini_time=0;//(b-b/ini);
        int N=v1;
        int a=N;
        int b=nc.GetNumRecs();
        int total_size=v2;
        vector <double> MSD (total_size, 0);
	vector <double> MSD_x (total_size, 0);
	vector <double> MSD_y (total_size, 0);
        //vector <double> Q (total_size, 0);
        double tau_samples=0;
        vector <int> t0 (total_size, 0);
        int size_MSD=b*N;
        //vector <double> MSD_i (size_MSD,0);
        int amostras=v17;
        //int start=1;
        int fim=start;
        double tau_a=-1.0;
        vector <double> time(b,0);
        double max=v6;
	double min=max-v11;
        int time_delta=0;
        string filename;
        double dt=v5;
	int LA=last*amostras;
	//	vector <double> MSD_x_2 (LA*total_size, 0);
	//	vector <double> t0_2 (LA*total_size, 0);

	Dscalar boxsize=v14;
	Dscalar CellRadius=v15;

        EOMPtr spp = make_shared<selfPropelledParticleDynamics>(N);
        shared_ptr<VoronoiQuadraticEnergyNotConfluent> spv  = make_shared<VoronoiQuadraticEnergyNotConfluent>(N, false, boxsize);
        spv->setCellRadiusUniform(v15);
        spv->setCellAdaptationTimeUniform(v10);
        SimulationPtr sim = make_shared<Simulation>();
        sim->setConfiguration(spv);
        sim->addUpdater(spp,spv);
        bool initializeGPU=false;
        sim->setCPUOperation(!initializeGPU);
        FILE *f;
 
        for(int avg=(start-1)*amostras+1; avg<=fim*amostras; avg++){
            cont=cont+1;
            vector <double> unrappedposx(N,0);
            vector <double> unrappedposy(N,0);
            
            //vector <vector <double> > vardata(2*N, vector<double> (b, 0));
            //double vardata[2*N][b];
            double L;
            //vector <vector <double> > vararea(2*N, vector<double> (b, 0));
            //double vararea[2*N][b];
            //vector <vector <double> > varper(N, vector<double> (b, 0));
            //double varper[N][b];
            //double val1=0;
            //double val2=0;
            //double val3=0;
            //double val4=0;

            int curr_line=(avg-1)/amostras+1;           
            //int v1, v2, v3, v4, v9, v12, v14;
            //double v5, v6, v7, v8, v10, v11, v13;
            //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/dados.txt";
	    //string access_dados="/media/hdd/Results/scripts25/dados.txt";
	    //string access_dados="/media/diogo/Elements/Tissues_results/scripts/dados2.txt";
	    file.open(access_dados);

            for(int line=0; line<curr_line; line++){
                file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11 >> v12 >> v13 >> v14 >> v15 >> v16 >> v17;
            }

            file.close();

            string add_JOB_num=to_string(avg);
            //if(avg<10)add_JOB_num="00"+to_string(avg);
            //else if(avg<100)add_JOB_num="0"+to_string(avg);
            //else add_JOB_num=to_string(avg);
            
            cout<<add_JOB_num<<endl;
           
            char dataname[1000000]; 
            //sprintf(dataname, "/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/scripts39/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/GPU/scripts5/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
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
                ncdat.ReadState(spv,0,false);
                b=ncdat.GetNumRecs();
                //printf("b=%u\n", b);
            }
           
            ArrayHandle<Dscalar2> h_p(spv->cellPositions,access_location::host,access_mode::read);
            ArrayHandle<Dscalar2> h_ap(spv->AreaPeriPreferences,access_location::host,access_mode::read);
            ArrayHandle<Dscalar2> h_AP(spv->AreaPeri,access_location::host,access_mode::read);

            Dscalar x11,x12,x21,x22;
            spv->Box->getBoxDims(x11,x12,x21,x22);
            L=x11;

            /*for(int r=0; r<b; r++){
	      if(r>0)ncdat.ReadState(spv,r,false);
	      time[r]= spv->currentTime; 
	      int rrr=0;
	      for (int rr = 0; rr < N; ++rr){
                
		vardata[rrr][r]=h_p.data[rr].x;
		vardata[rrr+1][r]=h_p.data[rr].y;
                
		//vararea[rr][r]=h_ap.data[rr].x;
		//vararea[N+rr][r]=h_AP.data[rr].x;
                    
		varper[rr][r]=h_ap.data[rr].y;
                    
		rrr=rrr+2;
	      }
	      }*/
	    
  
            //int vec[6]={1350, 1500, 1650, 1700, 1850, 2000};
            //int ini_time=(b-b/ini);
            //int ini_time=vec[ini-1];
	    	    
            if((cont-1)%amostras==0 || cont==1){
	      fill(MSD.begin(), MSD.end(), 0);
	      fill(MSD_x.begin(), MSD_x.end(), 0);
	      fill(MSD_y.begin(), MSD_y.end(), 0);	
	      //fill(Q.begin(), Q.end(), 0);
	      fill(t0.begin(), t0.end(), 0);
	      //fill(MSD_i.begin(), MSD_i.end(), 0);
	      t0[0]=1.0;
	      tau_a=-1;
            }

            //Qtest=zeros([b-ini_time 1]);

            int tau_bool=-1;
	    int start_nn=1;
	    ncdat.ReadState(spv,start_nn,false);
	    time[start_nn]= spv->currentTime;
            vector <double> xlast(N, 0);
            vector <double> ylast(N, 0);
	    vector <double> time_d(N, time[start_nn]);
	    vector <int> time_n(N, 0);
	    
            for(int nn = start_nn; nn<start_nn+1; nn++){

	      if(nn>start_nn)
		{
		  ncdat.ReadState(spv,nn,false);
		  time[nn]= spv->currentTime;
		}
	      
                vector <double> xn(N, 0);
                vector <double> yn(N, 0);
		vector <int> type(N, 0);
		
                double xm = 0.0;
                double ym = 0.0;
                int site=0;
                for(int cc=0; cc<N; cc=cc+1){
		  xn[site] = h_p.data[cc].x;
		  yn[site] = h_p.data[cc].y;
                    
		  xlast[site]=xn[site];
		  ylast[site]=yn[site];
                    
		  unrappedposy[site]=yn[site];
		  unrappedposx[site]=xn[site];
                    
		  xm=xm+xn[site];
		  ym=ym+yn[site];

		  /*if(h_ap.data[cc].y>max-0.00001)*/type[cc]=1;
                    
		  site++;
                }
                xm=xm/N;
                ym=ym/N;
    
                /*double xm = 0.0;
                double ym = 0.0;
                double max_dim=N;
        
                for(int ii=0; ii<max_dim; ii++){
                    xm=xm+xn[ii];
                    ym=ym+yn[ii];
                }
                xm=xm/max_dim;
                ym=ym/max_dim;*/
		
                for(int n = nn+1; n<b; n++){
                    
		  //time_delta=(time[n]-time[nn])/dt;
		  
                    ncdat.ReadState(spv,n,false);
		    time[n]= spv->currentTime;
		    
                    vector <double> x(N, 0);
                    vector <double> y(N, 0);
                    double xmi=0;
                    double ymi=0;
                    site=0;
		    double save_site=0;
                    for(int cc=0; cc<N; cc=cc+1){
		      save_site=x[site];
		      x[site] = h_p.data[cc].x;
		      y[site] = h_p.data[cc].y;
			
		      xmi=xmi+x[site];
		      ymi=ymi+y[site];
                        
		      site++;
                    }
                    xmi=xmi/N;
                    ymi=ymi/N;
                    
                    
                    /*double xmi=0;
                    double ymi=0;
            
                    for (int ii=0; ii<N; ii++){
                        xmi=xmi+x[ii];
                        ymi=ymi+y[ii];
                    }
                    xmi=xmi/N;
                    ymi=ymi/N;*/
                    
                    double posxm=abs(xmi-xm);
                    double posym=abs(ymi-ym);
    
                    if(posxm>L/2){
                        posxm=L-posxm;
                    }
                    if(posym>L/2){
                        posym=L-posym;
                    }
    
                    double maxVal=-100000;
                    int diff_N=0;
                    
                    int MSD_final=0;
                    for(int ii=0; ii<N; ii++){

		      /* if(varper[ii][n]>3.9-0.00001 && type[ii]==0)
			{
			  type[ii]=1;
			  xn[ii]=x[ii];
			  yn[ii]=y[ii];
			  xm=xmi;
			  ym=ymi;
			  time_d[ii]=time[n];
			  unrappedposx[ii]=x[ii];
			  unrappedposy[ii]=y[ii];
			  time_n[ii]=n;
			  }*/
		      
		      //if(h_ap.data[ii].y<max-0.00001)type[ii]=0;
		      
		      if(type[ii]==1){
			  
			time_delta=(time[n]-time_d[ii])/dt;
			  
			double posx=x[ii]-xlast[ii];
			double posy=y[ii]-ylast[ii];
			if(posx>L/2)posx=L-posx;
			if(posx<-L/2)posx=-L-posx;
			if(posy>L/2)posy=L-posy;
			if(posy<-L/2)posy=-L-posy;
			unrappedposx[ii]=unrappedposx[ii]+posx;
			unrappedposy[ii]=unrappedposy[ii]+posy;
			posx=unrappedposx[ii]-xn[ii];
			posy=unrappedposy[ii]-yn[ii];
			posx=posx-posxm;
			posy=posy-posym;
			double final_calc=posx*posx+posy*posy;
			diff_N++;
                            
			int idx=ii+(n-time_n[ii])*N;
			//MSD_i[idx]=MSD_i[idx]+final_calc;
			MSD[time_delta]=MSD[time_delta]+final_calc;
			MSD_x[time_delta]=MSD_x[time_delta]+(posx*posx);
			MSD_y[time_delta]=MSD_y[time_delta]+(posy*posy);
			t0[time_delta]=t0[time_delta]+1;

			//MSD_x_2[(avg-1)+time_delta*LA]+=(posx*posx);
			//t0_2[(avg-1)+time_delta*LA]++;
		      
			MSD_final+=final_calc;
			if(maxVal<final_calc/N)maxVal=final_calc;
			//if(sqrt(final_calc)<sqrt(vararea[ii][n])/2)Q[time_delta]=Q[time_delta]+1;
		      }	

		      xlast[ii]=x[ii];
		      ylast[ii]=y[ii];	
                    }

		    //MSD_x_s[avg]=MSD_x_2

		    /*if(diff_N>0)
		      {
			MSD[time_delta]=MSD[time_delta]+MSD_final/diff_N;
			t0[time_delta]=t0[time_delta]+1;
			}*/
		    
		    //else MSD[time_delta]=0;
		    
                    //if(diff_N>0){
                      //  MSD[time_delta]=MSD[time_delta]/diff_N;
                      //  Q[time_delta]=Q[time_delta]/diff_N;
                    //}
                }
            }

            if(cont%amostras==0 && cont>0){
                int value=cont/amostras;
		//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Stripe/s39/MSD_"+to_string(start) +".txt";
		//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Stripe/GPU/s5/MSD_"+to_string(start)+".txt";
		filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Non_Confluent/scripts"+to_string(scriptsIni)+"/MSD_"+to_string(start)+".txt";
                ofstream file1;
                file1.open (filename);
                for(int n = 1; n<total_size; n++){

		  if(MSD_x[n]>0 /*|| MSD_y[n]>0*/){
		      //t0[n]=t0[n];
                        //Q[n]=Q[n]/N;
                        //time[n]=time[n]-(10000*0.001);
            
                        MSD[n]=MSD[n]/t0[n];
			MSD_x[n]=MSD_x[n]/t0[n];
			
			double MSD_2=0;
			for(int kk=avg-amostras; kk<avg; kk++)
			  {
			    //MSD_2+=((MSD_x_2[kk+n*LA]/t0_2[kk+n*LA])-MSD_x[n])*((MSD_x_2[kk+n*LA]/t0_2[kk+n*LA])-MSD_x[n]);
			  }
			MSD_2/=(amostras-1);
			
			MSD_y[n]=MSD_y[n]/t0[n];
                        //Q[n]=Q[n]/t0[n];
                    
                        //if(Q[n]<exp(-1) && tau_a<0 && n>0)tau_a=n*dt;
			
                        file1 << MSD[n] <<" "<< sqrt(MSD_2)/sqrt(amostras)  <<" "<< MSD_x[n] <<" "<<  MSD_y[n] << " " << "Q[n]" << " " << n*dt << " " << tau_a << " "<< t0[n] <<endl;
                    }
                }
                file1.close();

		/*
                //filename="/home/destevao/CellGPUSubstrate/Results/"+to_string(ini)+"/data_distribution_"+to_string(value)+".txt";
                //ofstream file1
		filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/data_MSD_"+to_string(start) +".txt";
                file1.open (filename);
                for(int n = 0; n<b*N; n++){
                    int y=n/N;
                    int x=n%N;
                    //if(MSD_i[n]){
                
		    //file1 << MSD_i[n]/amostras << " " << time[y]/0.001 << " " << x << " " << y << endl;

                    //}
                }
                file1.close();*/
            }
        }
        
        
    }
    return 0;
}
