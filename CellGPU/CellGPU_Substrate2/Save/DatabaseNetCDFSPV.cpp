#include "DatabaseNetCDFSPV.h"
/*! \file DatabaseNetCDFSPV.cpp */

SPVDatabaseNetCDF::SPVDatabaseNetCDF(int np, string fn, NcFile::FileMode mode, bool exclude)
    : BaseDatabaseNetCDF(fn,mode),
      Nv(np),
      Current(0),
      exclusions(exclude)
    {
    switch(Mode)
        {
        case NcFile::ReadOnly:
            break;
        case NcFile::Write:
            GetDimVar();
            break;
        case NcFile::Replace:
            SetDimVar();
            break;
        case NcFile::New:
            SetDimVar();
            break;
        default:
            ;
        };
    }

void SPVDatabaseNetCDF::SetDimVar()
    {
    //Set the dimensions
    recDim = File.add_dim("rec");
    NvDim  = File.add_dim("Nv",  Nv);
    dofDim = File.add_dim("dof", Nv*2);
    boxDim = File.add_dim("boxdim",4);
    unitDim = File.add_dim("unit",1);

    //Set the variables
    timeVar          = File.add_var("time",     ncDscalar,recDim, unitDim);
    means0Var          = File.add_var("means0",     ncDscalar,recDim, unitDim);
    forceMaxVar          = File.add_var("ForceMax",     ncDscalar,recDim, unitDim);
    posVar          = File.add_var("pos",       ncDscalar,recDim, dofDim);
    typeVar          = File.add_var("type",         ncInt,recDim, NvDim );
    directorVar          = File.add_var("director",         ncDscalar,recDim, NvDim );
    BoxMatrixVar    = File.add_var("BoxMatrix", ncDscalar,recDim, boxDim);
    PrefPerimeterVar    = File.add_var("PrefPerimeter", ncDscalar,recDim, NvDim);
    PrefAreaVar    = File.add_var("PrefArea", ncDscalar,recDim, dofDim);
    PrefSelfPropVar    = File.add_var("PrefSelfProp", ncDscalar,recDim, NvDim);
    if(exclusions)
        exVar          = File.add_var("externalForce",       ncDscalar,recDim, dofDim);
    }

void SPVDatabaseNetCDF::GetDimVar()
    {
    //Get the dimensions
    recDim = File.get_dim("rec");
    boxDim = File.get_dim("boxdim");
    NvDim  = File.get_dim("Nv");
    dofDim = File.get_dim("dof");
    unitDim = File.get_dim("unit");
    //Get the variables
    posVar          = File.get_var("pos");
    typeVar          = File.get_var("type");
    directorVar          = File.get_var("director");
    means0Var          = File.get_var("means0");
    forceMaxVar          = File.get_var("ForceMax");
    BoxMatrixVar    = File.get_var("BoxMatrix");
    timeVar    = File.get_var("time");
    PrefPerimeterVar    = File.get_var("PrefPerimeter");
    PrefAreaVar    = File.get_var("PrefArea");
    PrefSelfPropVar    = File.get_var("PrefSelfProp");
    if(exclusions)
        exVar = File.get_var("externalForce");
    }

void SPVDatabaseNetCDF::WriteState(STATE s, Dscalar time, int rec)
    {
    if(rec<0)   rec = recDim->size();
    if (time < 0) time = s->currentTime;

    std::vector<Dscalar> boxdat(4,0.0);
    Dscalar x11,x12,x21,x22;
    s->Box->getBoxDims(x11,x12,x21,x22);
    boxdat[0]=x11;
    boxdat[1]=x12;
    boxdat[2]=x21;
    boxdat[3]=x22;

    std::vector<Dscalar> posdat(2*Nv);
    std::vector<Dscalar> directordat(Nv);
    std::vector<int> typedat(Nv);
    std::vector<Dscalar> p0(Nv);
    std::vector<Dscalar> a0(2*Nv);
    std::vector<Dscalar> v0(Nv);
    int idx = 0;
    Dscalar means0=0.0;

    Dscalar fmax=0.0;
    ArrayHandle<Dscalar2> h_f(s->cellForces,access_location::host,access_mode::read);
    Dscalar fx = 0.0;
    Dscalar fy = 0.0;
    Dscalar min = 10000;
    Dscalar max = -10000;   
 
    ArrayHandle<Dscalar2> h_p(s->cellPositions,access_location::host,access_mode::read);
    ArrayHandle<Dscalar> h_cd(s->cellDirectors,access_location::host,access_mode::read);
    ArrayHandle<int> h_ct(s->cellType,access_location::host,access_mode::read);
    ArrayHandle<int> h_ex(s->exclusions,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_ap(s->AreaPeriPreferences,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_mot(s->Motility,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_AP(s->AreaPeri,access_location::host,access_mode::read);

    for (int ii = 0; ii < Nv; ++ii)
        {
        int pidx = s->tagToIdx[ii];
        Dscalar px = h_p.data[pidx].x;
        Dscalar py = h_p.data[pidx].y;
        posdat[(2*idx)] = px;
        posdat[(2*idx)+1] = py;
        p0[ii]=h_ap.data[pidx].y;
        a0[ii]=h_ap.data[pidx].x;
        a0[Nv+ii]=h_AP.data[pidx].x;
        v0[ii]=h_mot.data[pidx].x;
        directordat[ii] = h_cd.data[pidx];
        if(h_ex.data[ii] == 0)
            typedat[ii] = h_ct.data[pidx];
        else
            typedat[ii] = h_ct.data[pidx]-5;
        idx +=1;

        if (h_f.data[pidx].y >max)
            max = h_f.data[pidx].y;
        if (h_f.data[pidx].x >max)
            max = h_f.data[pidx].x;
        if (h_f.data[pidx].y < min)
            min = h_f.data[pidx].y;
        if (h_f.data[pidx].x < min)
            min = h_f.data[pidx].x;
        fx += h_f.data[pidx].x;
        fy += h_f.data[pidx].y;
        };
    if(max>abs(min))fmax=max;
    else fmax=abs(min);

//    means0 = means0/Nv;
    means0 = s->reportq();

    //Write all the data
    means0Var      ->put_rec(&means0,      rec);
    forceMaxVar      ->put_rec(&fmax,      rec);
    timeVar      ->put_rec(&time,      rec);
    posVar      ->put_rec(&posdat[0],     rec);
    typeVar       ->put_rec(&typedat[0],      rec);
    directorVar       ->put_rec(&directordat[0],      rec);
    BoxMatrixVar->put_rec(&boxdat[0],     rec);
    PrefPerimeterVar  ->  put_rec(&p0[0],   rec);    
    PrefAreaVar  ->  put_rec(&a0[0],   rec);
    PrefSelfPropVar  ->  put_rec(&v0[0],   rec);

    if(exclusions)
        {
        ArrayHandle<Dscalar2> h_ef(s->external_forces,access_location::host,access_mode::read);
        std::vector<Dscalar> exdat(2*Nv);
        int id = 0;
        for (int ii = 0; ii < Nv; ++ii)
            {
            int pidx = s->tagToIdx[ii];
            Dscalar px = h_ef.data[pidx].x;
            Dscalar py = h_ef.data[pidx].y;
            exdat[(2*id)] = px;
            exdat[(2*id)+1] = py;
            id +=1;
            };
        exVar      ->put_rec(&exdat[0],     rec);
        };

    File.sync();
    }

void SPVDatabaseNetCDF::ReadState(STATE t, int rec,bool geometry)
    {

    //initialize the NetCDF dimensions and variables
    //test if there is exclusion data to read...
    int tester = File.num_vars();
    if (tester == 7)
        exclusions = true;
    GetDimVar();

    //get the current time
    timeVar-> set_cur(rec);
    timeVar->get(& t->currentTime,1,1);

    //set the box
    BoxMatrixVar-> set_cur(rec);
    std::vector<Dscalar> boxdata(4,0.0);
    BoxMatrixVar->get(&boxdata[0],1, boxDim->size());
    t->Box->setGeneral(boxdata[0],boxdata[1],boxdata[2],boxdata[3]);

    //get the positions
    posVar-> set_cur(rec);
    std::vector<Dscalar> posdata(2*Nv,0.0);
    posVar->get(&posdata[0],1, dofDim->size());

    ArrayHandle<Dscalar2> h_p(t->cellPositions,access_location::host,access_mode::overwrite);
    for (int idx = 0; idx < Nv; ++idx)
        {
        Dscalar px = posdata[(2*idx)];
        Dscalar py = posdata[(2*idx)+1];
        h_p.data[idx].x=px;
        h_p.data[idx].y=py;
        };

    //get the properties of the cells
    PrefPerimeterVar  ->  set_cur(rec);
    PrefAreaVar  ->  set_cur(rec);
    PrefSelfPropVar  ->  set_cur(rec);
    std::vector<Dscalar> p0data(Nv, 0.0);
    std::vector<Dscalar> a0data(2*Nv, 0.0);
    std::vector<Dscalar> v0data(Nv, 0.0);
    PrefPerimeterVar  ->  get(&p0data[0],1,   NvDim->size());
    PrefAreaVar  ->  get(&a0data[0],1,   dofDim->size());
    PrefSelfPropVar  ->  get(&v0data[0],1,   NvDim->size());

    ArrayHandle<Dscalar2> h_ap(t->AreaPeriPreferences,access_location::host,access_mode::overwrite);
    ArrayHandle<Dscalar2> h_mot(t->Motility,access_location::host,access_mode::overwrite);
    ArrayHandle<Dscalar2> h_AP(t->AreaPeri,access_location::host,access_mode::overwrite);

    for (int idx = 0; idx < Nv; ++idx){
        Dscalar a0val = a0data[idx];
        Dscalar p0val = p0data[idx];
        Dscalar aval = a0data[Nv+idx];
        Dscalar v0val = v0data[idx];

        h_ap.data[idx].x=a0val;
        h_ap.data[idx].y=p0val;
        h_AP.data[idx].x=aval;
        h_mot.data[idx].x=v0val;
    };


    //get cell types and cell directors
    typeVar->set_cur(rec);
    std::vector<int> ctdata(Nv,0.0);
    typeVar->get(&ctdata[0],1, NvDim->size());
    ArrayHandle<int> h_ct(t->cellType,access_location::host,access_mode::overwrite);

    directorVar->set_cur(rec);
    std::vector<Dscalar> cddata(Nv,0.0);
    directorVar->get(&cddata[0],1, NvDim->size());
    ArrayHandle<Dscalar> h_cd(t->cellDirectors,access_location::host,access_mode::overwrite);
    for (int idx = 0; idx < Nv; ++idx)
        {
        h_cd.data[idx]=cddata[idx];;
        h_ct.data[idx]=ctdata[idx];;
        };

    //read in excluded forces if applicable...
    if (tester == 7)
        {
        exVar-> set_cur(rec);
        std::vector<Dscalar> efdata(2*Nv,0.0);
        exVar->get(&posdata[0],1, dofDim->size());
        ArrayHandle<Dscalar2> h_ef(t->external_forces,access_location::host,access_mode::overwrite);
        for (int idx = 0; idx < Nv; ++idx)
            {
            Dscalar efx = efdata[(2*idx)];
            Dscalar efy = efdata[(2*idx)+1];
            h_ef.data[idx].x=efx;
            h_ef.data[idx].y=efy;
            };
        };


    //by default, compute the triangulation and geometrical information
    /*
    if(geometry)
        {
        t->resetDelLocPoints();
        t->updateCellList();
        t->globalTriangulationCGAL();
        t->resetLists();
        if(t->GPUcompute)
            t->computeGeometryGPU();
        else
            t->computeGeometryCPU();
        };
    */
    }


