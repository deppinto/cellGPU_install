
#include "voronoiQuadraticEnergyNotConfluent.h"
#include "cuda_profiler_api.h"
/*! \file voronoiQuadraticEnergyNotConfluent.cpp */

/*!
\param n number of cells to initialize
\param reprod should the simulation be reproducible (i.e. call a RNG with a fixed seed)
\post initializeVoronoiQuadraticEnergy(n,initGPURNcellsG) is called, as is setCellPreferenceUniform(1.0,4.0)
*/

VoronoiQuadraticEnergyNotConfluent::VoronoiQuadraticEnergyNotConfluent(int n, bool reprod, Dscalar boxsize)
    {
    printf("Initializing %i cells with random positions in a square box... \n",n);
    Reproducible = reprod;
    setBoxSize(boxsize);
    initializeVoronoiQuadraticEnergyNotConfluent(n);
    setCellPreferencesUniform(1.0,4.0);
    };

/*!
\param  n Number of cells to initialized
\post all GPUArrays are set to the correct size, v0 is set to 0.05, Dr is set to 1.0, the
Hilbert sorting period is set to -1 (i.e. off), the moduli are set to KA=KP=1.0, voronoiModelBase is
initialized (initializeVoronoiModelBase(n) gets called), particle exclusions are turned off, and auxiliary
data structures for the topology are set
*/
//take care of all class initialization functions
void VoronoiQuadraticEnergyNotConfluent::initializeVoronoiQuadraticEnergyNotConfluent(int n)
    {
    initializeVoronoiModelBase(n);
    Timestep = 0;
    setDeltaT(0.01);
    };

/*!
goes through the process of computing the forces on either the CPU or GPU, either with or without
exclusions, as determined by the flags. Assumes the geometry has NOT yet been computed.
\post the geometry is computed, and force per cell is computed.
*/
void VoronoiQuadraticEnergyNotConfluent::computeForces()
    {
    if(forcesUpToDate)
       return; 
    forcesUpToDate = true;
    computeGeometry();  
    for (int ii = 0; ii < Ncells; ++ii)
      computeVoronoiNonConfluentForceCPU(ii);

    //exit (911);
    };


void VoronoiQuadraticEnergyNotConfluent::reportError(ofstream &outfile)
    {
    ArrayHandle<Dscalar2> p(cellPositions,access_location::host,access_mode::read);
    //outfile << Ncells <<endl;
    for (int ii = 0; ii < Ncells ; ++ii)
        outfile << p.data[ii].x <<"\t" <<p.data[ii].y <<endl;
    };


/*!
\param i The particle index for which to compute the net force, assuming addition tension terms between unlike particles
\post the net force on cell i is computed
*/
void VoronoiQuadraticEnergyNotConfluent::computeVoronoiNonConfluentForceCPU(int i)
    {
      int qq=-1;
      int tstep=-1;

    //read in all the data we'll need
    ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_f(cellForces,access_location::host,access_mode::readwrite);
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_APpref(AreaPeriPreferences,access_location::host,access_mode::read);
    ArrayHandle<Dscalar> h_l(cellRadius,access_location::host,access_mode::read);
	
    ArrayHandle<int> h_ht(voroType,access_location::host,access_mode::read);
    ArrayHandle<int> h_vn(voroNum,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_v(NCvoroCur,access_location::host,access_mode::read);
    ArrayHandle<Dscalar4> h_t(voroTriangulation,access_location::host,access_mode::read);
    ArrayHandle<int> h_vcn(voroCellNeigh,access_location::host,access_mode::read);
    ArrayHandle<Dscalar> h_teta(voroOutAngle,access_location::host,access_mode::read);

    ArrayHandle<Dscalar2> h_external_forces(external_forces,access_location::host,access_mode::overwrite);
    ArrayHandle<int> h_exes(exclusions,access_location::host,access_mode::read);

    nc_idx = Index2D(2*neighMax,Ncells);
    
    if(i==qq && timestep==tstep)cout<<"START:------------------ "<<qq<<" "<<timestep<<endl;

    //get Delaunay neighbors of the cell
    int neigh = h_vn.data[i];
    vector<Dscalar4> ns(neigh);
    vector<int> h_type(neigh);
    vector<Dscalar2> voro(neigh);
    vector<Dscalar> angle(neigh);
    for (int nn = 0; nn < neigh; ++nn)
      {
	ns[nn]=h_t.data[nc_idx(nn,i)];
	h_type[nn]=h_ht.data[nc_idx(nn,i)];
	voro[nn]=h_v.data[nc_idx(nn,i)];
	angle[nn]=h_teta.data[nc_idx(nn,i)];
      };
    
    //compute base set of voronoi points, and the derivatives of those points w/r/t cell i's position
    vector<Matrix2x2> dhdri(neigh);
    Matrix2x2 Id;
    Dscalar2 rij,rik, rjk;
    Dscalar2 pi = h_p.data[i];
    Dscalar l=h_l.data[i];
    Dscalar2 forceSum;
    forceSum.x=0.0;forceSum.y=0.0;
    Dscalar Adiff = KA*(h_AP.data[i].x - h_APpref.data[i].x);
    Dscalar Pdiff = KP*(h_AP.data[i].y - h_APpref.data[i].y);
    Dscalar2 vlast,vnext,vother, vcur;

    for (int nn = 0; nn < neigh;++nn)
      {

	if(h_type[nn]==1)
	  {
	    rij.x=ns[nn].x;
	    rij.y=ns[nn].y;
	    rik.x=ns[nn].z;
	    rik.y=ns[nn].w;
	    rjk.x =rik.x-rij.x;
	    rjk.y =rik.y-rij.y;
	    
	    Dscalar betaD = -dot(rik,rik)*dot(rij,rjk);
	    Dscalar gammaD = dot(rij,rij)*dot(rik,rjk);
	    Dscalar alphaD = dot(rjk, rjk)*dot(rij, rik);
	    Dscalar cp = rij.x*rjk.y - rij.y*rjk.x;
	    Dscalar D = 2*cp*cp;
	    Dscalar2 dbDdri,dgDdri,dDdriOD,z;

	    z.x = betaD*rij.x+gammaD*rik.x;
	    z.y = betaD*rij.y+gammaD*rik.y;

	    dbDdri.x = 2*dot(rij,rjk)*rik.x+dot(rik,rik)*rjk.x;
	    dbDdri.y = 2*dot(rij,rjk)*rik.y+dot(rik,rik)*rjk.y;

	    dgDdri.x = -2*dot(rik,rjk)*rij.x-dot(rij,rij)*rjk.x;
	    dgDdri.y = -2*dot(rik,rjk)*rij.y-dot(rij,rij)*rjk.y;

	    dDdriOD.x = (-2.0*rjk.y)/cp;
	    dDdriOD.y = (2.0*rjk.x)/cp;

	    dhdri[nn] = Id+1.0/D*(dyad(rij,dbDdri)+dyad(rik,dgDdri)-(betaD+gammaD)*Id-dyad(z,dDdriOD));

	    if(i==qq && timestep==tstep)cout<<"here1: "<<dhdri[nn].x11<<" "<<dhdri[nn].x12<<" "<<dhdri[nn].x21<<" "<<dhdri[nn].x22<<endl;
	  }
	else
	  {
	    rij.x=ns[nn].x;
	    rij.y=ns[nn].y;
	    int signt=ns[nn].z;
	    Dscalar lambda=dot(rij,rij);

	    Dscalar Dval=2*sqrt(4*l*l-lambda)*lambda*sqrt(lambda);
	    dhdri[nn].x11=0.5+signt*(-4*l*l*rij.x*rij.y)/Dval;
	    dhdri[nn].x21=signt*(-4*l*l*lambda+lambda*lambda+(rij.x)*(rij.x)*4*l*l)/Dval;
	    dhdri[nn].x12=signt*(4*l*l*lambda-lambda*lambda-(rij.y)*(rij.y)*4*l*l)/Dval;
	    dhdri[nn].x22=0.5+signt*(4*l*l*rij.x*rij.y)/Dval;
	    if(i==qq && timestep==tstep)cout<<"HERE: "<< signt <<" "<< sqrt(rij.x*rij.x+rij.y*rij.y) <<" "<<rij.x<<" "<<rij.y<<" "<<lambda<<" "<<Dval<<" "<<dhdri[nn].x11<<" "<<dhdri[nn].x12<<" "<<dhdri[nn].x21<<" "<<dhdri[nn].x22 <<endl;
	  }

	//start calculating forces
//if(i==qq)l=5.0;
	int save_prev_idx=nn-1;
        if(save_prev_idx<0)save_prev_idx+=neigh;
	vlast = voro[save_prev_idx];

	//first, let's do the self-term, dE_i/dr_i
	vcur = voro[nn];
	int val=(nn+1)%neigh;
	vnext = voro[val];

	Dscalar2 dAidv,dPidv,dAm,dAn,dPm,dPn;
	bool signArea;
	AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[val], vcur, vnext, l, false, angle[nn]);
	AreaPeriCalc_for_Forces(dPn, dAn, h_type[save_prev_idx], h_type[nn], vcur, vlast, l, true, angle[save_prev_idx]);
	dAidv = dAm + dAn;
	dPidv = dPm + dPn;

	if(i==qq && timestep==tstep)cout<<"here2: "<<dPidv.x<<" "<<dPidv.y<<" "<<dAidv.x<<" "<<dAidv.y<<" "<<Adiff<<" "<<Pdiff<<endl;
	if(i==qq && timestep==tstep)cout<<setprecision(15)<<"here0: "<<vcur.x<<" "<<vcur.y<<" "<<vnext.x<<" "<<vnext.y<<" "<<vlast.x<<" "<<vlast.y<<" "<<h_p.data[i].x<<" "<<h_p.data[i].y<<endl;
	if(i==qq && timestep==tstep)cout<<"here00 "<<  h_type[nn]<<" "<<h_type[val] <<" "<< h_type[save_prev_idx]<<" "<<angle[nn]<<" "<<angle[save_prev_idx]<<" "<<((vnext.x*vcur.x+vnext.y*vcur.y)*(vnext.x*vcur.x+vnext.y*vcur.y)/(l*l*l*l))<<" "<<sin(angle[nn]) <<" "<<sin(angle[save_prev_idx])<<endl;
	//now let's compute the other terms...first we need to find the third voronoi
	//position that v_cur is connected to

	int baseNeigh=h_vcn.data[nc_idx(nn, i)];
        int lastNeigh=h_vcn.data[nc_idx(save_prev_idx, i)];
	int neigh2=h_vn.data[baseNeigh];
	int DT_other_idx=-1;
	int nnType = h_type[nn];
	int otherNeigh;
	Dscalar2 rjh;
	int savePoint=h_vcn.data[nc_idx(neigh2-1,baseNeigh)];

	Dscalar2 dAjdv;
	dAjdv.x=0;
	dAjdv.y=0;
	Dscalar2 dPjdv;
	dPjdv.x=0;
	dPjdv.y=0;
	Dscalar Pjdiff=0;
	Dscalar Ajdiff=0;
	Dscalar2 dAkdv,dPkdv;

	Dscalar Akdiff = KA*(h_AP.data[baseNeigh].x  - h_APpref.data[baseNeigh].x);
	Dscalar Pkdiff = KP*(h_AP.data[baseNeigh].y  - h_APpref.data[baseNeigh].y);
	int testPoint;
	int testType;

	if(nnType>1)
	  {
	    for (int n2 = 0; n2 < neigh2; ++n2)
	      {
		testPoint = h_vcn.data[nc_idx(n2,baseNeigh)];
		testType = h_ht.data[nc_idx(n2,baseNeigh)];
		switch(nnType)
		  {
		  case 2:
		    if(testPoint == i && testType==2)
		      {
                        if(savePoint==i)
                          {
			    DT_other_idx = (n2+1)%neigh2;
		  	    int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
			    Dscalar otherAngle=h_teta.data[nc_idx(n2,baseNeigh)];
			    int prev_other=n2-1;
			    if(prev_other<0)prev_other+=neigh2;
			    Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			    vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			    rjh=h_v.data[nc_idx(n2,baseNeigh)];

		  	    AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], otherType, rjh, vother, l, false, otherAngle);	    
		  	    //AreaPeriCalc_for_Forces(dPn, dAn, h_type[val], h_type[nn], vcur, vnext, l, true, 0.0);
			    AreaPeriCalc_for_Forces(dPn, dAn, h_type[val], h_type[nn], rjh, jnext, l, true, 0.0);
		  	    dAkdv = dAm + dAn;
		  	    dPkdv = dPm + dPn;
			    n2=neigh2;
                          }
                        else
                          {
			    DT_other_idx = n2-1;
                            if(DT_other_idx<0)DT_other_idx=neigh2+DT_other_idx;
		  	    int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
			    Dscalar otherAngle=h_teta.data[nc_idx(DT_other_idx,baseNeigh)];
			    int prev_other=(n2+1)%neigh2;
			    Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			    vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			    rjh=h_v.data[nc_idx(n2,baseNeigh)];

	    	  	    //AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[save_prev_idx], vcur, vlast, l, false, 0.0);	 
	    	  	    AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[save_prev_idx], rjh, jnext, l, false, 0.0);   
	    	  	    AreaPeriCalc_for_Forces(dPn, dAn, otherType, h_type[nn], rjh, vother, l, true, otherAngle);
		  	    dAkdv = dAm + dAn;
		  	    dPkdv = dPm + dPn;
			    n2=neigh2;
                          }
		      }
		    break;
		  case 3:
		    if(testPoint == i && testType==4)
		      {
			DT_other_idx = (n2+1)%neigh2;
		  	int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
		        Dscalar otherAngle=h_teta.data[nc_idx(n2,baseNeigh)];
			int prev_other=n2-1;
			if(prev_other<0)prev_other+=neigh2;
			Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			rjh=h_v.data[nc_idx(n2,baseNeigh)];

		  	AreaPeriCalc_for_Forces(dPm, dAm, 4, otherType, rjh, vother, l, false, otherAngle);	    
		  	//AreaPeriCalc_for_Forces(dPn, dAn, 3, 4, vcur, vnext, l, true, 0.0);
			AreaPeriCalc_for_Forces(dPn, dAn, 3, 4, rjh, jnext, l, true, 0.0);
		  	dAkdv = dAm + dAn;
		  	dPkdv = dPm + dPn;
			n2=neigh2;
		      }
		    break;
		  case 4:
		    if(testPoint == i && testType==3)
		      {
			DT_other_idx = (n2-1);
                        if(DT_other_idx<0)DT_other_idx=neigh2+DT_other_idx;
		  	int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
		        Dscalar otherAngle=h_teta.data[nc_idx(DT_other_idx,baseNeigh)];
			int prev_other=(n2+1)%neigh2;
			Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			rjh=h_v.data[nc_idx(n2,baseNeigh)];

	    	  	//AreaPeriCalc_for_Forces(dPm, dAm, 3, 4, vcur, vlast, l, false, 0.0);
			AreaPeriCalc_for_Forces(dPm, dAm, 3, 4, rjh, jnext, l, false, 0.0);	    
	    	  	AreaPeriCalc_for_Forces(dPn, dAn, otherType, 3, rjh, vother, l, true, otherAngle);
		  	dAkdv = dAm + dAn;
		  	dPkdv = dPm + dPn;
			n2=neigh2;
		      }
		    break;
		  }

                savePoint=testPoint;
	      };
	  }
	else
	  {
	    for (int n2 = 0; n2 < neigh2; ++n2)
	      {
		int testPoint = h_vcn.data[nc_idx(n2, baseNeigh)];
		if(testPoint == lastNeigh){
		  DT_other_idx = (n2+1)%neigh2;
		  //if(DT_other_idx<0)DT_other_idx=neigh2+DT_other_idx;
		  otherNeigh=h_vcn.data[nc_idx(DT_other_idx,baseNeigh)];
		  //DT_other_idx = DT_other_idx-1;
		  //if(DT_other_idx<0)DT_other_idx+=neigh2; 
		  rjh=h_v.data[nc_idx(DT_other_idx,baseNeigh)]+h_p.data[baseNeigh];
		  Box->putInBoxReal(rjh);

	          //Dscalar2 nn1 = h_p.data[baseNeigh];
		  //Dscalar2 r1;
		  Box->minDist(rjh,pi,vother);
		  //Box->minDist(rjh,r1,vother);
		  //vother=rjh+r1;
		  int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];

		  AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], otherType, vcur, vother, l, false, 0.0);	    
		  AreaPeriCalc_for_Forces(dPn, dAn, h_type[val], h_type[nn], vcur, vnext, l, true, 0.0);
		  dAkdv = dAm + dAn;
		  dPkdv = dPm + dPn;

		  Dscalar2 help;
		  Circumcenter(h_p.data[lastNeigh], h_p.data[baseNeigh], h_p.data[i], help);
		  if(i==qq && timestep==tstep)cout<<"FINAL DECIDER: "<<help.x<<" "<<help.y<<" "<<vcur.x<<" "<<vcur.y<<" "<<h_p.data[i].x<<" "<<h_p.data[i].y<<endl;
		  if(i==qq && timestep==tstep)cout<<"here01: "<<vother.x<<" "<<vother.y<<" "<</*r1.x<<" "<<r1.y<<" "<<*/rjh.x<<" "<<rjh.y<<" "<<n2<<" "<<neigh2<<endl;

	   	  Ajdiff = KA*(h_AP.data[lastNeigh].x - h_APpref.data[lastNeigh].x);
	    	  Pjdiff = KP*(h_AP.data[lastNeigh].y - h_APpref.data[lastNeigh].y);

	    	  AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[save_prev_idx], vcur, vlast, l, false, 0.0);	    
	    	  AreaPeriCalc_for_Forces(dPn, dAn, otherType, h_type[nn], vcur, vother, l, true, 0.0);
	    	  dAjdv = dAm + dAn;
	    	  dPjdv = dPm + dPn;

                  n2=neigh2;
		}
	      };
	  }
	
	if(DT_other_idx == -1)
	  {
	    printf("Triangulation problem %i\n",DT_other_idx);
	    cout<<i<<" "<<baseNeigh<<" "<<neigh<<" "<<neigh2<<" "<<nnType<<" "<<h_vcn.data[nc_idx((DT_other_idx+1)%neigh2,baseNeigh)] <<endl;
	    char fn[256];
	    sprintf(fn,"failed.txt");
	    ofstream output(fn);
	    reportError(output);
	    throw std::exception();
	  };
	

	if(i==qq && timestep==tstep)cout<<"here3: "<<Akdiff<<" "<<dAkdv.x<<" "<<Pkdiff<<" "<<dPkdv.x<<" "<<dAkdv.y<<" "<<dPkdv.y<<endl;
	if(i==qq && timestep==tstep)cout<<"here05: "<<baseNeigh<<" "<<otherNeigh<<" "<<lastNeigh<<endl;
	if(i==qq && timestep==tstep)cout<<"here06: "<<vother.x<<" "<<vother.y<<" "<< nnType<<" "<<testType<<" "<<testPoint <<endl;
	if(i==qq && timestep==tstep)cout<<"here4: "<<Ajdiff<<" "<<dAjdv.x<<" "<<Pjdiff<<" "<<dPjdv.x<<" "<<dAjdv.y<<" "<<dPjdv.y<<endl;

	//energy calculation
	Dscalar2 dEdv;

	dEdv.x = 2.0*Adiff*dAidv.x + 2.0*Pdiff*dPidv.x;
	dEdv.y = 2.0*Adiff*dAidv.y + 2.0*Pdiff*dPidv.y;
	if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Adiff*dAidv.x<<" "<<2.0*Pdiff*dPidv.x<<" "<<2.0*Adiff*dAidv.y<<" "<<2.0*Pdiff*dPidv.y<<endl;
	dEdv.x += 2.0*Akdiff*dAkdv.x + 2.0*Pkdiff*dPkdv.x;
	dEdv.y += 2.0*Akdiff*dAkdv.y + 2.0*Pkdiff*dPkdv.y;
	if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Akdiff*dAkdv.x<<" "<<2.0*Pkdiff*dPkdv.x<<" "<<2.0*Akdiff*dAkdv.y<<" "<<2.0*Pkdiff*dPkdv.y<<endl;
	dEdv.x += 2.0*Ajdiff*dAjdv.x + 2.0*Pjdiff*dPjdv.x;
	dEdv.y += 2.0*Ajdiff*dAjdv.y + 2.0*Pjdiff*dPjdv.y;
	if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Ajdiff*dAjdv.x<<" "<<2.0*Pjdiff*dPjdv.x<<" "<<2.0*Ajdiff*dAjdv.y<<" "<<2.0*Pjdiff*dPjdv.y<<endl;

	Dscalar2 temp = dEdv*dhdri[nn];
	forceSum.x += temp.x;
	forceSum.y += temp.y;

	if(i==qq && timestep==tstep)cout<<"FINAL: "<<temp.x<<" "<<temp.y<<" "<<dEdv.x<<" "<<dEdv.y<<" "<<forceSum.x<<" "<<forceSum.y <<endl;

	//take care of the explicit derivative of the energy only when curved segments exist
        bool curved=true;    
        if(h_type[save_prev_idx]==1 || h_type[nn]==1)curved=false;
        else if(h_type[save_prev_idx]==3 && h_type[nn]==4)curved=false;
	
	//The convention here is different so as to use the same functions as before
	if(curved==true)
	{
          Adiff = KA*(h_AP.data[i].x - h_APpref.data[i].x);
          Pdiff = KP*(h_AP.data[i].y - h_APpref.data[i].y);
	  AreaPeriCalc_for_Forces(dPm, dAm, h_type[save_prev_idx], h_type[nn], vcur, vlast, l, false, angle[save_prev_idx]);
	  AreaPeriCalc_for_Forces(dPn, dAn, h_type[save_prev_idx], h_type[nn], vlast, vcur, l, true, angle[save_prev_idx]);
	  dAidv = dAm + dAn;
	  dPidv = dPm + dPn;
	  if(i==qq && timestep==tstep)cout<<"AP1: "<<dAidv.x<<" "<<dPidv.x<<" "<<dAidv.y<<" "<<dPidv.y<<endl;
	  dEdv.x = 2.0*Adiff*dAidv.x + 2.0*Pdiff*dPidv.x;
	  dEdv.y = 2.0*Adiff*dAidv.y + 2.0*Pdiff*dPidv.y;
	  if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Adiff*dAidv.x<<" "<<2.0*Pdiff*dPidv.x<<" "<<2.0*Adiff*dAidv.y<<" "<<2.0*Pdiff*dPidv.y<<endl;
	  forceSum.x -= dEdv.x;
	  forceSum.y -= dEdv.y;
	  if(i==qq && timestep==tstep)cout<<"here002: "<<dPidv.x<<" "<<dPidv.y<<" "<<dAidv.x<<" "<<dAidv.y<<" "<<Adiff<<" "<<Pdiff<<endl;
	  if(i==qq && timestep==tstep)cout<<"FINAL2: "<<dEdv.x<<" "<<dEdv.y<<" "<<forceSum.x<<" "<<forceSum.y <<endl;
	}

	vlast=vcur;
	save_prev_idx=nn;
      };

    h_f.data[i].x=forceSum.x;
    h_f.data[i].y=forceSum.y;
    if(particleExclusions)
      {
        if(h_exes.data[i] != 0)
            {
            h_f.data[i].x = 0.0;
            h_f.data[i].y = 0.0;
            h_external_forces.data[i].x=-forceSum.x;
            h_external_forces.data[i].y=-forceSum.y;
            };
        }

    if(i==qq && timestep==tstep){
      cout<<"FORCES CALC NON CONFLUENT "<< i <<" "<<neigh<<" "<< h_f.data[i].x<<" "<< h_f.data[i].y <<" "<<h_AP.data[i].x <<" "<<h_AP.data[i].y <<endl;
    char fn[256];
    sprintf(fn,"snap001.pov");
    ofstream output(fn);
        if(timestep==tstep){writeTriangulation(output);exit (911);}
      }
    };

//the minus sign from the force calculations are introduced in this function by using the appropriate conventions
void VoronoiQuadraticEnergyNotConfluent::AreaPeriCalc_for_Forces(Dscalar2 &dPm, Dscalar2 &dAm, int type1, int type2, Dscalar2 vcur, Dscalar2 vnext, Dscalar l, bool sign, Dscalar teta)
    {

    Dscalar Pthreshold = THRESHOLD;
    bool curved=true;
    
    if(type1==1 || type2==1)curved=false;
    else if(type1==3 && type2==4)curved=false;

    //if(debug)cout<<"AREA_PERI_CALC: "<<curved<<endl;
    
    if(curved==false)
      {
	int signt=1;
	if(sign==false)signt=-1;

	dAm.x = 0.5*signt*vnext.y;
	dAm.y = 0.5*signt*(-vnext.x);

	Dscalar2 dnext;
	dnext.x = (vnext.x-vcur.x);
	dnext.y = (vnext.y-vcur.y);
	Dscalar dnnorm = sqrt(dnext.x*dnext.x+dnext.y*dnext.y);
	if(dnnorm < Pthreshold)
	  dnnorm = Pthreshold;

	dPm.x=dnext.x/dnnorm;
	dPm.y=dnext.y/dnnorm;

	//if(debug)cout<<"AREA_PERI_CALC: false; "<< dAm.x<<" "<<dAm.y<<" "<<dPm.x<<" "<<dPm.y  <<endl;
      }
    else
      {
	int signt=1;
	//if(teta>PI)signt=-1;

	Dscalar dnorm=l*l*sin(teta);
		//sqrt(1-((vnext.x*vcur.x+vnext.y*vcur.y)*(vnext.x*vcur.x+vnext.y*vcur.y)/(l*l*l*l)));	
        //if(debug)cout<<vcur.x<< " "<<vcur.y<<" "<<vnext.x<<" "<<vnext.y<<" "<<type1<<" "<<type2<<" "<<((vnext.x*vcur.x+vnext.y*vcur.y)*(vnext.x*vcur.x+vnext.y*vcur.y)/(l*l*l*l)) <<endl;
	dAm.x=signt*0.5*l*l*vnext.x/dnorm;
	dAm.y=signt*0.5*l*l*vnext.y/dnorm;

	dPm.x=signt*l*vnext.x/dnorm;
	dPm.y=signt*l*vnext.y/dnorm;

	//if(debug)cout<<"AREA_PERI_CALC: true; "<< dAm.x<<" "<<dAm.y<<" "<<dPm.x<<" "<<dPm.y<<" "<<teta  <<endl;
	//if(((vnext.x*vcur.x+vnext.y*vcur.y)*(vnext.x*vcur.x+vnext.y*vcur.y)/(l*l*l*l))>=1)exit(911);
      }
    }

/*!
Returns the quadratic energy functional:
E = \sum_{cells} K_A(A_i-A_i,0)^2 + K_P(P_i-P_i,0)^2
*/
Dscalar VoronoiQuadraticEnergyNotConfluent::computeEnergy()
    {
    if(!forcesUpToDate)
        computeForces();
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_APP(AreaPeriPreferences,access_location::host,access_mode::read);
    Energy = 0.0;
    for (int nn = 0; nn  < Ncells; ++nn)
        {
        Energy += KA * (h_AP.data[nn].x-h_APP.data[nn].x)*(h_AP.data[nn].x-h_APP.data[nn].x);
        Energy += KP * (h_AP.data[nn].y-h_APP.data[nn].y)*(h_AP.data[nn].y-h_APP.data[nn].y);
        };
    return Energy;
    };

void VoronoiQuadraticEnergyNotConfluent::computeVoronoiNonConfluentForceCPU(int i, int j, Dscalar2 &fij)
    {
      int qq=-1;
      int tstep=-1;

    //read in all the data we'll need
    ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_f(cellForces,access_location::host,access_mode::readwrite);
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_APpref(AreaPeriPreferences,access_location::host,access_mode::read);
    ArrayHandle<Dscalar> h_l(cellRadius,access_location::host,access_mode::read);
	
    ArrayHandle<int> h_ht(voroType,access_location::host,access_mode::read);
    ArrayHandle<int> h_vn(voroNum,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_v(NCvoroCur,access_location::host,access_mode::read);
    ArrayHandle<Dscalar4> h_t(voroTriangulation,access_location::host,access_mode::read);
    ArrayHandle<int> h_vcn(voroCellNeigh,access_location::host,access_mode::read);
    ArrayHandle<Dscalar> h_teta(voroOutAngle,access_location::host,access_mode::read);

    nc_idx = Index2D(2*neighMax,Ncells);
    
    if(i==qq && timestep==tstep)cout<<"START:------------------ "<<qq<<" "<<timestep<<endl;

    //get Delaunay neighbors of the cell
    int neigh = h_vn.data[i];
    vector<Dscalar4> ns(neigh);
    vector<int> h_type(neigh);
    vector<Dscalar2> voro(neigh);
    vector<Dscalar> angle(neigh);
    for (int nn = 0; nn < neigh; ++nn)
      {
	ns[nn]=h_t.data[nc_idx(nn,i)];
	h_type[nn]=h_ht.data[nc_idx(nn,i)];
	voro[nn]=h_v.data[nc_idx(nn,i)];
	angle[nn]=h_teta.data[nc_idx(nn,i)];
      };
    
    //compute base set of voronoi points, and the derivatives of those points w/r/t cell i's position
    vector<Matrix2x2> dhdri(neigh);
    Matrix2x2 Id;
    Dscalar2 rij,rik, rjk;
    Dscalar2 pi = h_p.data[i];
    Dscalar l=h_l.data[i];
    Dscalar2 forceSum;
    forceSum.x=0.0;forceSum.y=0.0;
    Dscalar Adiff = KA*(h_AP.data[i].x - h_APpref.data[i].x);
    Dscalar Pdiff = KP*(h_AP.data[i].y - h_APpref.data[i].y);
    Dscalar2 vlast,vnext,vother, vcur;

    for (int nn = j; nn < j+1;++nn)
      {
	if(h_type[nn]==1)
	  {
	    rij.x=ns[nn].x;
	    rij.y=ns[nn].y;
	    rik.x=ns[nn].z;
	    rik.y=ns[nn].w;
	    rjk.x =rik.x-rij.x;
	    rjk.y =rik.y-rij.y;
	    
	    Dscalar betaD = -dot(rik,rik)*dot(rij,rjk);
	    Dscalar gammaD = dot(rij,rij)*dot(rik,rjk);
	    Dscalar alphaD = dot(rjk, rjk)*dot(rij, rik);
	    Dscalar cp = rij.x*rjk.y - rij.y*rjk.x;
	    Dscalar D = 2*cp*cp;
	    Dscalar2 dbDdri,dgDdri,dDdriOD,z;

	    z.x = betaD*rij.x+gammaD*rik.x;
	    z.y = betaD*rij.y+gammaD*rik.y;

	    dbDdri.x = 2*dot(rij,rjk)*rik.x+dot(rik,rik)*rjk.x;
	    dbDdri.y = 2*dot(rij,rjk)*rik.y+dot(rik,rik)*rjk.y;

	    dgDdri.x = -2*dot(rik,rjk)*rij.x-dot(rij,rij)*rjk.x;
	    dgDdri.y = -2*dot(rik,rjk)*rij.y-dot(rij,rij)*rjk.y;

	    dDdriOD.x = (-2.0*rjk.y)/cp;
	    dDdriOD.y = (2.0*rjk.x)/cp;

	    dhdri[nn] = Id+1.0/D*(dyad(rij,dbDdri)+dyad(rik,dgDdri)-(betaD+gammaD)*Id-dyad(z,dDdriOD));

	    if(i==qq && timestep==tstep)cout<<"here1: "<<dhdri[nn].x11<<" "<<dhdri[nn].x12<<" "<<dhdri[nn].x21<<" "<<dhdri[nn].x22<<endl;
	  }
	else
	  {
	    rij.x=ns[nn].x;
	    rij.y=ns[nn].y;
	    int signt=ns[nn].z;
	    Dscalar lambda=dot(rij,rij);

	    Dscalar Dval=2*sqrt(4*l*l-lambda)*lambda*sqrt(lambda);
	    dhdri[nn].x11=0.5+signt*(-4*l*l*rij.x*rij.y)/Dval;
	    dhdri[nn].x21=signt*(-4*l*l*lambda+lambda*lambda+(rij.x)*(rij.x)*4*l*l)/Dval;
	    dhdri[nn].x12=signt*(4*l*l*lambda-lambda*lambda-(rij.y)*(rij.y)*4*l*l)/Dval;
	    dhdri[nn].x22=0.5+signt*(4*l*l*rij.x*rij.y)/Dval;
	    if(i==qq && timestep==tstep)cout<<"HERE: "<< signt <<" "<< sqrt(rij.x*rij.x+rij.y*rij.y) <<" "<<rij.x<<" "<<rij.y<<" "<<lambda<<" "<<Dval<<" "<<dhdri[nn].x11<<" "<<dhdri[nn].x12<<" "<<dhdri[nn].x21<<" "<<dhdri[nn].x22 <<endl;
	  }

	//start calculating forces
//if(i==qq)l=5.0;
	int save_prev_idx=nn-1;
        if(save_prev_idx<0)save_prev_idx+=neigh;
	vlast = voro[save_prev_idx];

	//first, let's do the self-term, dE_i/dr_i
	vcur = voro[nn];
	int val=(nn+1)%neigh;
	vnext = voro[val];

	Dscalar2 dAidv,dPidv,dAm,dAn,dPm,dPn;
	bool signArea;
	AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[val], vcur, vnext, l, false, angle[nn]);
	AreaPeriCalc_for_Forces(dPn, dAn, h_type[save_prev_idx], h_type[nn], vcur, vlast, l, true, angle[save_prev_idx]);
	dAidv = dAm + dAn;
	dPidv = dPm + dPn;

	if(i==qq && timestep==tstep)cout<<"here2: "<<dPidv.x<<" "<<dPidv.y<<" "<<dAidv.x<<" "<<dAidv.y<<" "<<Adiff<<" "<<Pdiff<<endl;
	if(i==qq && timestep==tstep)cout<<setprecision(15)<<"here0: "<<vcur.x<<" "<<vcur.y<<" "<<vnext.x<<" "<<vnext.y<<" "<<vlast.x<<" "<<vlast.y<<" "<<h_p.data[i].x<<" "<<h_p.data[i].y<<endl;
	if(i==qq && timestep==tstep)cout<<"here00 "<<  h_type[nn]<<" "<<h_type[val] <<" "<< h_type[save_prev_idx]<<" "<<angle[nn]<<" "<<angle[save_prev_idx]<<" "<<((vnext.x*vcur.x+vnext.y*vcur.y)*(vnext.x*vcur.x+vnext.y*vcur.y)/(l*l*l*l))<<" "<<sin(angle[nn]) <<" "<<sin(angle[save_prev_idx])<<endl;
	//now let's compute the other terms...first we need to find the third voronoi
	//position that v_cur is connected to

	int baseNeigh=h_vcn.data[nc_idx(nn, i)];
        int lastNeigh=h_vcn.data[nc_idx(save_prev_idx, i)];
	int neigh2=h_vn.data[baseNeigh];
	int DT_other_idx=-1;
	int nnType = h_type[nn];
	int otherNeigh;
	Dscalar2 rjh;
	int savePoint=h_vcn.data[nc_idx(neigh2-1,baseNeigh)];

	Dscalar2 dAjdv;
	dAjdv.x=0;
	dAjdv.y=0;
	Dscalar2 dPjdv;
	dPjdv.x=0;
	dPjdv.y=0;
	Dscalar Pjdiff=0;
	Dscalar Ajdiff=0;
	Dscalar2 dAkdv,dPkdv;

	Dscalar Akdiff = KA*(h_AP.data[baseNeigh].x  - h_APpref.data[baseNeigh].x);
	Dscalar Pkdiff = KP*(h_AP.data[baseNeigh].y  - h_APpref.data[baseNeigh].y);
	int testPoint;
	int testType;

	if(nnType>1)
	  {
	    for (int n2 = 0; n2 < neigh2; ++n2)
	      {
		testPoint = h_vcn.data[nc_idx(n2,baseNeigh)];
		testType = h_ht.data[nc_idx(n2,baseNeigh)];
		switch(nnType)
		  {
		  case 2:
		    if(testPoint == i && testType==2)
		      {
                        if(savePoint==i)
                          {
			    DT_other_idx = (n2+1)%neigh2;
		  	    int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
			    Dscalar otherAngle=h_teta.data[nc_idx(n2,baseNeigh)];
			    int prev_other=n2-1;
			    if(prev_other<0)prev_other+=neigh2;
			    Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			    vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			    rjh=h_v.data[nc_idx(n2,baseNeigh)];

		  	    AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], otherType, rjh, vother, l, false, otherAngle);	    
		  	    //AreaPeriCalc_for_Forces(dPn, dAn, h_type[val], h_type[nn], vcur, vnext, l, true, 0.0);
			    AreaPeriCalc_for_Forces(dPn, dAn, h_type[val], h_type[nn], rjh, jnext, l, true, 0.0);
		  	    dAkdv = dAm + dAn;
		  	    dPkdv = dPm + dPn;
			    n2=neigh2;
                          }
                        else
                          {
			    DT_other_idx = n2-1;
                            if(DT_other_idx<0)DT_other_idx=neigh2+DT_other_idx;
		  	    int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
			    Dscalar otherAngle=h_teta.data[nc_idx(DT_other_idx,baseNeigh)];
			    int prev_other=(n2+1)%neigh2;
			    Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			    vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			    rjh=h_v.data[nc_idx(n2,baseNeigh)];

	    	  	    //AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[save_prev_idx], vcur, vlast, l, false, 0.0);	 
	    	  	    AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[save_prev_idx], rjh, jnext, l, false, 0.0);   
	    	  	    AreaPeriCalc_for_Forces(dPn, dAn, otherType, h_type[nn], rjh, vother, l, true, otherAngle);
		  	    dAkdv = dAm + dAn;
		  	    dPkdv = dPm + dPn;
			    n2=neigh2;
                          }
		      }
		    break;
		  case 3:
		    if(testPoint == i && testType==4)
		      {
			DT_other_idx = (n2+1)%neigh2;
		  	int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
		        Dscalar otherAngle=h_teta.data[nc_idx(n2,baseNeigh)];
			int prev_other=n2-1;
			if(prev_other<0)prev_other+=neigh2;
			Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			rjh=h_v.data[nc_idx(n2,baseNeigh)];

		  	AreaPeriCalc_for_Forces(dPm, dAm, 4, otherType, rjh, vother, l, false, otherAngle);	    
		  	//AreaPeriCalc_for_Forces(dPn, dAn, 3, 4, vcur, vnext, l, true, 0.0);
			AreaPeriCalc_for_Forces(dPn, dAn, 3, 4, rjh, jnext, l, true, 0.0);
		  	dAkdv = dAm + dAn;
		  	dPkdv = dPm + dPn;
			n2=neigh2;
		      }
		    break;
		  case 4:
		    if(testPoint == i && testType==3)
		      {
			DT_other_idx = (n2-1);
                        if(DT_other_idx<0)DT_other_idx=neigh2+DT_other_idx;
		  	int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];
		        Dscalar otherAngle=h_teta.data[nc_idx(DT_other_idx,baseNeigh)];
			int prev_other=(n2+1)%neigh2;
			Dscalar2 jnext=h_v.data[nc_idx(prev_other,baseNeigh)];
			vother=h_v.data[nc_idx(DT_other_idx,baseNeigh)];
			rjh=h_v.data[nc_idx(n2,baseNeigh)];

	    	  	//AreaPeriCalc_for_Forces(dPm, dAm, 3, 4, vcur, vlast, l, false, 0.0);
			AreaPeriCalc_for_Forces(dPm, dAm, 3, 4, rjh, jnext, l, false, 0.0);	    
	    	  	AreaPeriCalc_for_Forces(dPn, dAn, otherType, 3, rjh, vother, l, true, otherAngle);
		  	dAkdv = dAm + dAn;
		  	dPkdv = dPm + dPn;
			n2=neigh2;
		      }
		    break;
		  }

                savePoint=testPoint;
	      };
	  }
	else
	  {
	    for (int n2 = 0; n2 < neigh2; ++n2)
	      {
		int testPoint = h_vcn.data[nc_idx(n2, baseNeigh)];
		if(testPoint == lastNeigh){
		  DT_other_idx = (n2+1)%neigh2;
		  //if(DT_other_idx<0)DT_other_idx=neigh2+DT_other_idx;
		  otherNeigh=h_vcn.data[nc_idx(DT_other_idx,baseNeigh)];
		  //DT_other_idx = DT_other_idx-1;
		  //if(DT_other_idx<0)DT_other_idx+=neigh2; 
		  rjh=h_v.data[nc_idx(DT_other_idx,baseNeigh)]+h_p.data[baseNeigh];
		  Box->putInBoxReal(rjh);

	          //Dscalar2 nn1 = h_p.data[baseNeigh];
		  //Dscalar2 r1;
		  Box->minDist(rjh,pi,vother);
		  //Box->minDist(rjh,r1,vother);
		  //vother=rjh+r1;
		  int otherType=h_ht.data[nc_idx(DT_other_idx,baseNeigh)];

		  AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], otherType, vcur, vother, l, false, 0.0);	    
		  AreaPeriCalc_for_Forces(dPn, dAn, h_type[val], h_type[nn], vcur, vnext, l, true, 0.0);
		  dAkdv = dAm + dAn;
		  dPkdv = dPm + dPn;

		  Dscalar2 help;
		  Circumcenter(h_p.data[lastNeigh], h_p.data[baseNeigh], h_p.data[i], help);
		  if(i==qq && timestep==tstep)cout<<"FINAL DECIDER: "<<help.x<<" "<<help.y<<" "<<vcur.x<<" "<<vcur.y<<" "<<h_p.data[i].x<<" "<<h_p.data[i].y<<endl;
		  if(i==qq && timestep==tstep)cout<<"here01: "<<vother.x<<" "<<vother.y<<" "<</*r1.x<<" "<<r1.y<<" "<<*/rjh.x<<" "<<rjh.y<<" "<<n2<<" "<<neigh2<<endl;

	   	  Ajdiff = KA*(h_AP.data[lastNeigh].x - h_APpref.data[lastNeigh].x);
	    	  Pjdiff = KP*(h_AP.data[lastNeigh].y - h_APpref.data[lastNeigh].y);

	    	  AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[save_prev_idx], vcur, vlast, l, false, 0.0);	    
	    	  AreaPeriCalc_for_Forces(dPn, dAn, otherType, h_type[nn], vcur, vother, l, true, 0.0);
	    	  dAjdv = dAm + dAn;
	    	  dPjdv = dPm + dPn;

                  n2=neigh2;
		}
	      };
	  }
	
	if(DT_other_idx == -1)
	  {
	    printf("Triangulation problem %i\n",DT_other_idx);
	    cout<<i<<" "<<baseNeigh<<" "<<neigh<<" "<<neigh2<<" "<<nnType<<" "<<h_vcn.data[nc_idx((DT_other_idx+1)%neigh2,baseNeigh)] <<endl;
	    char fn[256];
	    sprintf(fn,"failed.txt");
	    ofstream output(fn);
	    reportError(output);
	    throw std::exception();
	  };
	

	if(i==qq && timestep==tstep)cout<<"here3: "<<Akdiff<<" "<<dAkdv.x<<" "<<Pkdiff<<" "<<dPkdv.x<<" "<<dAkdv.y<<" "<<dPkdv.y<<endl;
	if(i==qq && timestep==tstep)cout<<"here05: "<<baseNeigh<<" "<<otherNeigh<<" "<<lastNeigh<<endl;
	if(i==qq && timestep==tstep)cout<<"here06: "<<vother.x<<" "<<vother.y<<" "<< nnType<<" "<<testType<<" "<<testPoint <<endl;
	if(i==qq && timestep==tstep)cout<<"here4: "<<Ajdiff<<" "<<dAjdv.x<<" "<<Pjdiff<<" "<<dPjdv.x<<" "<<dAjdv.y<<" "<<dPjdv.y<<endl;

	//energy calculation
	Dscalar2 dEdv;

	dEdv.x = 2.0*Adiff*dAidv.x + 2.0*Pdiff*dPidv.x;
	dEdv.y = 2.0*Adiff*dAidv.y + 2.0*Pdiff*dPidv.y;
	if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Adiff*dAidv.x<<" "<<2.0*Pdiff*dPidv.x<<" "<<2.0*Adiff*dAidv.y<<" "<<2.0*Pdiff*dPidv.y<<endl;
	dEdv.x += 2.0*Akdiff*dAkdv.x + 2.0*Pkdiff*dPkdv.x;
	dEdv.y += 2.0*Akdiff*dAkdv.y + 2.0*Pkdiff*dPkdv.y;
	if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Akdiff*dAkdv.x<<" "<<2.0*Pkdiff*dPkdv.x<<" "<<2.0*Akdiff*dAkdv.y<<" "<<2.0*Pkdiff*dPkdv.y<<endl;
	dEdv.x += 2.0*Ajdiff*dAjdv.x + 2.0*Pjdiff*dPjdv.x;
	dEdv.y += 2.0*Ajdiff*dAjdv.y + 2.0*Pjdiff*dPjdv.y;
	if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Ajdiff*dAjdv.x<<" "<<2.0*Pjdiff*dPjdv.x<<" "<<2.0*Ajdiff*dAjdv.y<<" "<<2.0*Pjdiff*dPjdv.y<<endl;

	Dscalar2 temp = dEdv*dhdri[nn];
	forceSum.x += temp.x;
	forceSum.y += temp.y;

	if(i==qq && timestep==tstep)cout<<"FINAL: "<<temp.x<<" "<<temp.y<<" "<<dEdv.x<<" "<<dEdv.y<<" "<<forceSum.x<<" "<<forceSum.y <<endl;

	//take care of the explicit derivative of the energy only when curved segments exist
        bool curved=true;    
        if(h_type[save_prev_idx]==1 || h_type[nn]==1)curved=false;
        else if(h_type[save_prev_idx]==3 && h_type[nn]==4)curved=false;
	
	//The convention here is different so as to use the same functions as before
	if(curved==true)
	{
          Adiff = KA*(h_AP.data[i].x - h_APpref.data[i].x);
          Pdiff = KP*(h_AP.data[i].y - h_APpref.data[i].y);
	  AreaPeriCalc_for_Forces(dPm, dAm, h_type[save_prev_idx], h_type[nn], vcur, vlast, l, false, angle[save_prev_idx]);
	  AreaPeriCalc_for_Forces(dPn, dAn, h_type[save_prev_idx], h_type[nn], vlast, vcur, l, true, angle[save_prev_idx]);
	  dAidv = dAm + dAn;
	  dPidv = dPm + dPn;
	  if(i==qq && timestep==tstep)cout<<"AP1: "<<dAidv.x<<" "<<dPidv.x<<" "<<dAidv.y<<" "<<dPidv.y<<endl;
	  dEdv.x = 2.0*Adiff*dAidv.x + 2.0*Pdiff*dPidv.x;
	  dEdv.y = 2.0*Adiff*dAidv.y + 2.0*Pdiff*dPidv.y;
	  if(i==qq && timestep==tstep)cout<<"Energy1: "<<2.0*Adiff*dAidv.x<<" "<<2.0*Pdiff*dPidv.x<<" "<<2.0*Adiff*dAidv.y<<" "<<2.0*Pdiff*dPidv.y<<endl;
	  forceSum.x -= dEdv.x;
	  forceSum.y -= dEdv.y;
	  if(i==qq && timestep==tstep)cout<<"here002: "<<dPidv.x<<" "<<dPidv.y<<" "<<dAidv.x<<" "<<dAidv.y<<" "<<Adiff<<" "<<Pdiff<<endl;
	  if(i==qq && timestep==tstep)cout<<"FINAL2: "<<dEdv.x<<" "<<dEdv.y<<" "<<forceSum.x<<" "<<forceSum.y <<endl;
	}

	vlast=vcur;
	save_prev_idx=nn;
      };

    fij.x=forceSum.x;
    fij.y=forceSum.y;

    if(i==qq && timestep==tstep){
      cout<<"FORCES CALC NON CONFLUENT "<< i <<" "<<neigh<<" "<< h_f.data[i].x<<" "<< h_f.data[i].y <<" "<<h_AP.data[i].x <<" "<<h_AP.data[i].y <<endl;
    char fn[256];
    sprintf(fn,"snap001.pov");
    ofstream output(fn);
        if(timestep==tstep){writeTriangulation(output);exit (911);}
      }
    };


void VoronoiQuadraticEnergyNotConfluent::computeLocalStressTensor(int i, Matrix2x2 &sint)
    {

    //read in all the data we'll need
    ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_f(cellForces,access_location::host,access_mode::readwrite);
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_APpref(AreaPeriPreferences,access_location::host,access_mode::read);
    ArrayHandle<Dscalar> h_l(cellRadius,access_location::host,access_mode::read);

    ArrayHandle<int> h_ht(voroType,access_location::host,access_mode::read);
    ArrayHandle<int> h_vn(voroNum,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_v(NCvoroCur,access_location::host,access_mode::read);
    ArrayHandle<Dscalar4> h_t(voroTriangulation,access_location::host,access_mode::read);
    ArrayHandle<int> h_vcn(voroCellNeigh,access_location::host,access_mode::read);
    ArrayHandle<Dscalar> h_teta(voroOutAngle,access_location::host,access_mode::read);

    nc_idx = Index2D(2*neighMax,Ncells);

    //get Delaunay neighbors of the cell
    Dscalar2 v=h_p.data[i];
    int neigh = h_vn.data[i];
    vector<Dscalar4> ns(neigh);
    vector<int> h_type(neigh);
    vector<Dscalar2> voro(neigh);
    vector<Dscalar> angle(neigh);
    for (int nn = 0; nn < neigh; ++nn)
      {
        ns[nn]=h_t.data[nc_idx(nn,i)];
        h_type[nn]=h_ht.data[nc_idx(nn,i)];
        voro[nn]=h_v.data[nc_idx(nn,i)];
        angle[nn]=h_teta.data[nc_idx(nn,i)];
      };

    Dscalar l=1.0;
    Matrix2x2 P;
    Matrix2x2 T;
    P=0;
    T=0;
    Dscalar Pressure=2*(h_AP.data[i].x-h_APpref.data[i].x);
    Dscalar Tension=2*(h_AP.data[i].y-h_APpref.data[i].y);

    for (int nn = 0; nn < neigh; ++nn)
    {
	      int save_prev_idx=nn-1;
              if(save_prev_idx<0)save_prev_idx+=neigh;
              Dscalar2 vlast = voro[save_prev_idx];
              Dscalar2 vcur = voro[nn];
              int val=(nn+1)%neigh;
              Dscalar2 vnext = voro[val];

              Dscalar2 dAidv,dPidv,dAm,dAn,dPm,dPn;
              AreaPeriCalc_for_Forces(dPm, dAm, h_type[nn], h_type[val], vcur, vnext, l, false, angle[nn]);
              AreaPeriCalc_for_Forces(dPn, dAn, h_type[save_prev_idx], h_type[nn], vcur, vlast, l, true, angle[save_prev_idx]);
              dAidv = dAm + dAn;
              dPidv = dPm + dPn;

	      P+=dyad(voro[nn], dAidv);
	      T+=dyad(voro[nn], dPidv);
    }

    sint=Pressure*P+Tension*T;
    sint*=(1/h_AP.data[i].x);
    }

void VoronoiQuadraticEnergyNotConfluent::computeSwimStressTensor(int i, Matrix2x2 &sswim)
    {

    //read in all the data we'll need
    ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_APpref(AreaPeriPreferences,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_motility(Motility,access_location::host,access_mode::read);

    ArrayHandle<Dscalar> h_cd(cellDirectors, access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_v(cellVelocities, access_location::host,access_mode::read);

    Dscalar2 ni;
    ni.x=cos(h_cd.data[i]);
    ni.y=sin(h_cd.data[i]);
    Dscalar val=h_motility.data[i].x/h_AP.data[i].x;
    sswim.x11=val*ni.x*h_p.data[i].x;
    sswim.x12=val*ni.x*h_p.data[i].y;
    sswim.x21=val*ni.y*h_p.data[i].x;
    sswim.x22=val*ni.y*h_p.data[i].y;
    }
