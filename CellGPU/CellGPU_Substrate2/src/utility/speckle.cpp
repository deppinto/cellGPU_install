
#include "speckle.h"

////////////////////////////////////////////////////////////////////////
////								      //
////		    speckle_algorithm_fft  			      //
////								      //
////////////////////////////////////////////////////////////////////////


void speckle::fft(int sign, vector<complex<double> > &zs){
    unsigned int j = 0;
    int sizeZ = zs.size();
    
    for(unsigned int i = 0; i < sizeZ - 1; ++i){
        if(i < j) swap(zs[i], zs[j]);
        
        int m = sizeZ / 2;
        j ^= m;
        while((j & m) == 0){m /= 2; j ^= m;}
    }
    
    for(unsigned int j = 1; j < sizeZ; j *= 2){
        for(unsigned int m = 0;  m < j; ++m) {
            double t = M_PI * sign * m / j;
            complex<double> w = complex<double>(cos(t), sin(t));
            
            for(unsigned int i = m; i < sizeZ; i += 2 * j) {
                complex<double> zi = zs[i], t = w * zs.at(i + j);
                zs[i] = zi + t;
                zs.at(i + j) = zi - t;
            }
        }
    }
}


////////////////////////////////////////////////////////////////////////
////								      //
////		    speckle_algorithm_normal  			      //
////								      //
////////////////////////////////////////////////////////////////////////

double speckle::normal(double factor){
 //   noise.Reproducible = Reproducible;
    double rr1 = noise.getRealUniform(1e-6, 1.0);
    double rr2 = noise.getRealUniform(1e-6, 1.0);
    return sqrt(factor * log(rr1)) * cos(PI2 * rr2);
}


////////////////////////////////////////////////////////////////////////
////								      //
////		    speckle_random_field  			      //
////								      //
////////////////////////////////////////////////////////////////////////

void speckle::randomField(vector<vector<double> > &randField, double L, int nDisc, double sigma){
    
    //noise.Reproducible = Reproducible;
    double sigma2 = sigma * sigma;
    int N = nDisc;
    double delta_x = L / (1. * N);
    double factor = -2.;
    
    vector<vector<vector<double> > > kapas(N, vector<vector<double> >(N, vector<double>(2,0)));
    vector<vector<double> > as(N, vector<double>(N));
    vector<vector<complex<double> > > zs(N, vector<complex<double> >(N));
    
    double delta_k = 1. / L;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            kapas[i][j][0] = delta_k * (i - .5 * N);
            kapas[i][j][1] = delta_k * (j - .5 * N);
        }
    }
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            as[i][j] = sqrt(exp(- M_PI*M_PI*(kapas[i][j][0] * kapas[i][j][0] + kapas[i][j][1] * kapas[i][j][1]) * (2. * sigma2)));
        }
    }
    
    for(int i = 0; i <= .5 * N; i++){
        for(int j = 0; j <= .5 * N; j++){
            double R = as[i][j] * normal(factor);
            double phi = PI2 * noise.getRealUniform(1e-6, 1.0);
            //zs[i][j] = R * cos(phi) + 1i * R * sin(phi);
            zs[i][j] = complex<double>(R * cos(phi), R * sin(phi));
            
            int si, sj;
            if(i == 0){si = 0;}else{si = N - i;}
            if(j == 0){sj = 0;}else{sj = N - j;}
            zs[si][sj] = conj(zs[i][j]);
        }
    }
    
    zs[.5 * N][0]      = real(zs[.5 * N][0]);
    zs[0][.5 * N]      = real(zs[0][.5 * N]);
    zs[.5 * N][.5 * N] = real(zs[.5 * N][.5 * N]);
    
    for(int i = 1; i < .5 * N; i++){
        for(int j = 1; j < .5 * N; j++){
            double R = as[i][j] * normal(factor);
            double phi = PI2 * noise.getRealUniform(1e-6, 1.0);
            
            int si = N - i;
            int sj = N - j;
            //zs[i][sj] = R * cos(phi) + 1i * R * sin(phi);
            zs[i][sj] = complex<double>(R * cos(phi), R * sin(phi));
            zs[si][j] = conj(zs[i][sj]);
        }
    }

    for(int i = 0; i < N; i++){
        fft(-1, zs[i]);
    }
    vector<vector<complex<double> > > zs2(N, vector<complex<double> >(N));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            zs2[j][i] = zs[i][j];
        }
    }
    for(int j = 0; j < N; j++){
        fft(-1, zs2[j]);
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            zs[i][j] = zs2[j][i];
        }
    }
    
    double maxi = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            randField[i][j] = real(zs[i][j]);
            if((i % 2 == 0) && (j % 2 == 0)) randField[i][j] *= -1;
            if((i % 2 == 1) && (j % 2 == 1)) randField[i][j] *= -1;
            if(randField[i][j] > maxi) maxi = randField[i][j];
        }
    }
    
    for(int i = 0; i < 0.5 * N; i++){
        for(int j = 0; j < N; j++){
            int si = 0.5 * N + i;
            swap(randField[i][j],randField[si][j]);
        }
    }
    for(int i = 0; i < 0.5 * N; i++){
        for(int j = 0; j < N; j++){
            int si = .5 * N + i;
            swap(randField[j][i],randField[j][si]);
        }
    }
    
   /* 
    vector<double> list1(N*N);
    for(int i = 0; i < N * N; i++){
        list1[i] = -log(noise.getRealUniform(1e-6, 1.0));
    }
    sort(list1.begin(), list1.end());
    
    set< pair<double, int> > listExp;
    set< pair<double, int> >::iterator itListExp;
    
    int ii = 0;
    for(int v = 0; v < N; v++){
        for(int u = 0; u < N; u++){
            listExp.insert(make_pair(randField[v][u], ii));
            ii++;
        }
    }
    
    ii = 0;
    itListExp = listExp.begin();
    while(itListExp != listExp.end()){
        int pos = (*itListExp).second;
        randField[pos / N][pos % N] = list1[ii];
        ii++;
        itListExp++;
    }
    */

    maxi = 0;
    double mini=0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(randField[i][j] > maxi) maxi = randField[i][j];
            if(randField[i][j] < mini) mini = randField[i][j];
        }
    }
    
    
    ofstream file1;
    file1.open("speckle.txt");

    mini=abs(mini);
    double Delta=maxi+mini;

    for(int v = 0; v < N; v++){
        for(int u = 0; u < N; u++){
            randField[v][u] = randField[v][u] / Delta;
            file1 << randField[v][u] << " ";
        }
        file1 << endl;
    }
    file1.close();

   /*
 *      ofstream file_in1;
 *           file_in1.open("speckle.bin");
 *                for(int i = 0; i < N; i++){
 *                     file_in1.write(reinterpret_cast<char*>(&randField[i][0]), N * sizeof(double));
 *                          }
 *                               file_in1.close ();
 *                                    */

}

