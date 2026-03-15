#include "stl.h"
#include "count.h"
#include "sampling.h"
using namespace std;

// "sreg" is one subregion in a partition, linxi
void Count(shrink_reg_data_short & sreg, int* Tmp_count_train, int* Tmp_count_vali, bool* Tmp_unique, int dim, bool discrete, parameters & para){
    double* Tmp_upbound;
    Tmp_upbound =new double[dim];
    //int Tmp_count[dim];
    usint_mask Tmp_reg_code;
    //bool Tmp_unique[dim];
    for(int d=0; d< dim; d++){
        Tmp_count_train[d]=0;
		Tmp_count_vali[d]=0;
        Tmp_reg_code.x = (sreg.reg_code[d].x << 1); // this is the left/bottom piece, linxi
        Tmp_reg_code.mask = (sreg.reg_code[d].mask << 1) + 1;
        pair<double, double> ranges = convert_ranges(Tmp_reg_code);
        Tmp_upbound[d]=ranges.second;
        if(discrete) Tmp_unique[d] = true;
        else         Tmp_unique[d] = false;
    }
    for(int i=0; i<sreg.num; i++){
        for(int j=0; j< dim; j++){
            //if(sreg.regdata[i][j] <= Tmp_upbound[j]){
            if(para.data1D[para.pt_start[sreg.region_start+i]+j] <= Tmp_upbound[j]){
                ++Tmp_count_train[j];
            }
            //if(Tmp_unique[j] && fabs(sreg.regdata[0][j] - sreg.regdata[i][j]) > 1e-10)
            if(Tmp_unique[j] && fabs(para.data1D[para.pt_start[sreg.region_start]+j] - para.data1D[para.pt_start[sreg.region_start+i]+j]) > 1e-10)
                Tmp_unique[j] = (Tmp_unique[j] && false);
        }
    } 
    for(int i=0; i<sreg.valiNum; i++){
        for(int j=0; j< dim; j++){
            //if(sreg.regdata[i][j] <= Tmp_upbound[j]){
            if(para.data1D_vali[para.pt_start_vali[sreg.region_start_vali+i]+j] <= Tmp_upbound[j]){
                ++Tmp_count_vali[j];
            }
            //if(Tmp_unique[j] && fabs(sreg.regdata[0][j] - sreg.regdata[i][j]) > 1e-10)
            if(Tmp_unique[j] && fabs(para.data1D_vali[para.pt_start_vali[sreg.region_start_vali]+j] - para.data1D_vali[para.pt_start_vali[sreg.region_start_vali+i]+j]) > 1e-10)
                Tmp_unique[j] = (Tmp_unique[j] && false);
        }
    } 
    delete [] Tmp_upbound;
	/*
	for(int j=0; j <dim; j++) {
		cout<<Tmp_count_train[j]<<"\t"<<Tmp_count_vali[j]<<endl;
	}
	*/
}

// when calculating the weight, a normamlized version is used to simplify calculation, linxi
// this is BSP proposal
void Calculate_Weight_Gibbs(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count_train, int trainsize, int* Tmp_count_vali, int valisize, int RegionID){
    for(int i=0; i < dim; ++i){
        c_original[dim*RegionID+i] = 0;
        double beta=0.5; // a pseudo count, to avoid 0 count. 
        if ((double) trainsize / 200.0 < 0.5) beta = max(0.1, (double) trainsize / 200.0); // heuristic pseudo ct.
        if(Tmp_unique[i])
            c_original[dim*RegionID+i] = -MAXNUM;
        else{
            c_original[dim*RegionID+i] += ((double)Tmp_count_vali[i]) * log(beta+(double)Tmp_count_train[i]);
            c_original[dim*RegionID+i] += ((double)(valisize-Tmp_count_vali[i])) * log(beta+(double)(trainsize -Tmp_count_train[i]));
        }
        if(c_original[dim*RegionID+i] > 1 - MAXNUM){
            c_original[dim*RegionID+i] -= ((double)valisize) * log(2*beta+(double)trainsize);
            c_original[dim*RegionID+i] += (double)valisize * log(2.0);
        }
		//cout<<c_original[dim*RegionID+i]<<endl;
    }
}

// when calculating the weight, a normamlized version is used to simplify calculation, linxi
// this is BSP proposal
void Calculate_Weight_OPT(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count_train, int trainsize, int* Tmp_count_vali, int valisize, int RegionID){
    for(int i=0; i < dim; ++i){
        c_original[dim*RegionID+i] = 0;
        double beta=0.5; // a pseudo count, to avoid 0 count. 
        if ((double) trainsize / 200.0 < 0.5) beta = max(0.1, (double) trainsize / 200.0); // heuristic pseudo ct.
        if(Tmp_unique[i])
            c_original[dim*RegionID+i] = -MAXNUM;
        else{
            c_original[dim*RegionID+i] += ((double)Tmp_count_train[i]) * log(beta+(double)Tmp_count_train[i]);
            c_original[dim*RegionID+i] += ((double)(trainsize-Tmp_count_train[i])) * log(beta+(double)(trainsize -Tmp_count_train[i]));
        }
        if(c_original[dim*RegionID+i] > 1 - MAXNUM){
            c_original[dim*RegionID+i] -= ((double)trainsize) * log(2*beta+(double)trainsize);
            c_original[dim*RegionID+i] += (double)trainsize * log(2.0);
        }
		//cout<<c_original[dim*RegionID+i]<<endl;
    }
}

// when calculating the weight, a normamlized version is used to simplify calculation, linxi
// this is BSP proposal
void Calculate_Weight_RF(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count_train, int trainsize, int* Tmp_count_vali, int valisize, int RegionID){
    for(int i=0; i < dim; ++i){
		// uniform weight, randomly choose a direction to split. 
        c_original[dim*RegionID+i] = 0;
    }
}