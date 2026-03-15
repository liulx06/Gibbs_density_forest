

#include "examples.h"


// compare with john's method, most simple discontinuous case
void PK(){      //don't been called anywhere, thchiu

    parameters p;
    p.dim=3;
    p.levels=300;
    p.maxpercentage = 0.9;
    p.n = 100;
    p.resampling=2;
    p.samplesize=50000;
    p.steps=4;

    string temp = "result/PK.txt";
    freopen(temp.c_str(), "w", stdout);
    cout<<"START"<<endl;
    vector< map< OnePartition , double, CompairSRegs> > pathctmap;
    readdata("data/dim3_2norm_50000.txt", p.samplesize, p.dim, p.data, true);
    vector<OnePartition_data> Ps;
    int bestind;
    SIS(p, false, Ps, bestind);

     print_partition(Ps[bestind]);
     fclose(stdout);

}

void cytometry(){   //don't been called anywhere, thchiu

    parameters p;
    p.dim=15;
    p.levels=1000;
    p.maxpercentage = 0.9;
    p.n = 100;
    p.resampling=2;
    p.samplesize=212250;
    p.steps=4;

    string temp = "result/cyto" +  toStr<int>(p.levels) + "full.txt";
    freopen(temp.c_str(), "w", stdout);
    cout<<"START"<<endl;
    vector< map< OnePartition , double, CompairSRegs> > pathctmap;
    readdata("data/CyTOF54_Tube01_Day1_Unstim1_curated.fcs_eventnum_Ungated_Jnstim1_Day1_normalized_1_Unstim1_Singlets.txt", p.samplesize, p.dim, p.data, true);
    FSIS2( p,false, NULL,false, false, 0.001, pathctmap);

     fclose(stdout);

}


void cash(){    //don't been called anywhere, thchiu

    parameters p;
 //   p.dim=15;
    p.levels=100;
    p.maxpercentage = 0.9;
    p.n = 100;
    p.resampling=2;
    p.samplesize=3791;
    p.steps=4;
    vector<vector<double> > origdata;
    string temp = "result/0314_00006_" +  toStr<int>(p.levels) + ".txt";
    freopen(temp.c_str(), "w", stdout);
    cout<<"START"<<endl;
    vector< map< OnePartition , double, CompairSRegs> > pathctmap;
    readdata("data/0314_00006.txt", p.samplesize, 40, origdata, true);
    cout<<origdata[0][0]<<" "<<origdata[0][1]<<" "<<origdata[0][2]<<" "<<origdata[0][3]<<endl;
    vector<int> useind(17,0);
    useind[0] = 0; useind[1] = 1; useind[2] = 2; useind[3] = 5; useind[4] = 6;
    useind[5] = 7; useind[6] = 8; useind[7] = 9; useind[8] = 10; useind[9] = 11;
    useind[10] = 12; useind[11] = 13; useind[12] = 24; useind[13] = 27; useind[14] = 30;
    useind[15] = 33; useind[16] = 36;
    for(int i=0; i< p.samplesize; i++){
        vector<double> onepiece((int)useind.size(), 0);
        for(int j=0; j<(int)useind.size(); j++) {
            onepiece[j]=origdata[i][useind[j]+2];
 //           cout<<origdata[i][useind[j]+2]<<" ";
        }
//        cout<<"\n";
        p.data.push_back(onepiece);
    }
    p.dim = (int)useind.size();
    cout<<"dim="<<p.dim<<" "<<p.data[0][0]<<" "<<p.data[0][1]<<" "<<p.data[0][2]<<" "<<p.data[0][3]<<endl;
    FSIS2( p,false, NULL,false, true, 0, pathctmap);
    fclose(stdout);

}


//  density estimation problem, using copula. but here the density is the one after copula
void FSIS2(parameters & p,bool knowndensity, double (*density)(const vector<double> &x),   bool  pathctcorrect, bool discrete, double smoothneighbordist, vector< map< OnePartition , double, CompairSRegs> >& pathctmap){

    int dim = p.dim;
    int samplesize = p.samplesize;
    parameters para;
    para.maxpercentage = 0.9;
    para.samplesize = samplesize;
    para.steps = 5;
    para.dim = 1;
    para.n = 100;
    para.levels = 20;
    para.resampling = 2;

    vector<OnePartition_data> Ps;
    int bestind;
    vector<vector<double> >  dataF= p.data;
    for(int d=0; d<dim; d++){
        cout<<"dim="<<d<<endl<<endl;
        para.data.clear();
        vector<double> onepiece(1,0);
        vector<double> onevariable(samplesize, 0);
        for(int i=0; i<samplesize; i++)   {
            onepiece[0]=p.data[i][d];
            onevariable[i]= p.data[i][d];
            para.data.push_back(onepiece);
        }
        Ps.clear();
        SIS(p, discrete, Ps, bestind);
        OnePartition_data P= Ps[bestind];
        print_partition(P);
        vector<double> newx = Ftransform(onevariable, P);
        for(int i=0; i<samplesize; i++)        dataF[i][d]=newx[i];
    }
     para.dim = dim;
     para.data.clear();
     para.data=dataF;
     para.levels = p.levels;
     para.n = p.n;
     para.steps = p.steps;
     para.resampling = p.resampling;
     Ps.clear();
     SIS(p, discrete, Ps, bestind);

     Output_distance (  p,Ps,density, smoothneighbordist,pathctmap, pathctcorrect);

     print_partition(Ps[bestind]);
     return;

}


