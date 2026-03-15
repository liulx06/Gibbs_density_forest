#ifndef COUNT_H
#define COUNT_H
#include "sampling.h"
void Count(shrink_reg_data_short & sreg, int* Tmp_count_train, int* Tmp_count_vali, bool* Tmp_unique, int dim, bool discrete, parameters & p);
void Calculate_Weight_Gibbs(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count_train, int trainsize, int* Tmp_count_vali, int valisize, int RegionID);
void Calculate_Weight_OPT(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count_train, int trainsize, int* Tmp_count_vali, int valisize, int RegionID);
void Calculate_Weight_RF(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count_train, int trainsize, int* Tmp_count_vali, int valisize, int RegionID);
#endif
