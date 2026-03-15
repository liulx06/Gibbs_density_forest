#include "randomshift.h"
#include "timer_ed520.h"


    vector<double> marginalmax(k, -100000000);
    vector<double> marginalmin(k,  1000000000);
    for(int i=0; i<k; i++){
        for (int j=0; j<(int)data.size(); j++){
            if(data[j][i]>marginalmax[i])  marginalmax[i]=data[j][i];
            if(data[j][i]<marginalmin[i])  marginalmin[i]=data[j][i];

        }
        cerr<<"max="<<marginalmax[i]<<"; min="<<marginalmin[i]<<endl;
    }
    int sumsmall=0;
    for (int i = 0; i < k; i++) {
        if (marginalmin[i] == marginalmax[i]) {
            cout << "only one value in dim=" << i << endl;
            for (int j = 0; j < (int) data.size(); j++) {
                 data[j][i] = 0.45;
            }    
            continue;
        }
        for (int j = 0; j < (int) data.size(); j++) {
  //          data[j][i] = (data[j][i] - marginalmin[i] + 0.001) / (marginalmax[i] - marginalmin[i] + 0.002);
            data[j][i] = (data[j][i] - marginalmin[i]) / (marginalmax[i] - marginalmin[i]);
            if(data[j][i]<0) data[j][i]=0;
            if(data[j][i]>1) data[j][i]=1;

//            cerr<<data[j][i]<<'\t';
            if(data[j][i]<0.01) sumsmall+=1;

        }
    }
    cerr<<"sumsmall="<<sumsmall<<endl;
    mmax = marginalmax;
    mmin = marginalmin;