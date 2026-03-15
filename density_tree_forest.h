#ifndef DENSITY_TREE_FOREST_H
#define	DENSITY_TREE_FOREST_H

//#include "hellingerdist.h"
#include "sampling.h"
#include "SISfunctions.h"
#include "readdata.h"
#include "timer_ed520.h"    //+ by ed520

// function to fit Bayesian density forests under shifts, using Bayesian model averaging
void bma_density_forest(parameters & p, parameters & para, vector<OneDensity_data>& ensemble, int ensemble_sizeP, int n_shiftM, double shift_sizeM, int n_shiftP, double shift_sizeP, double betaM, double betaP, string ofilename);

// function to draw samples from the Gibbs posterior
void gibbs_density_forest(parameters & p, parameters & para, vector<OneForest_data>& forest_of_Fs, int depth_max, int depth_min, int ensemble_sizeP, int number_of_forests, double betaM, double betaP, double splitP, string ofilename, double tree_temperature = 1.0, double forest_temperature = 1.0, int Cdepth=30, bool UQ_flag=false);
// int n_shiftM, double shift_sizeM, int n_shiftP, double shift_sizeP,

// function to draw samples under the optional polya tree prior, with a copula transformation first
void optional_polya_tree_copula(parameters & p, parameters & para, vector<OnePartition_data>& marginalP, OnePartition_data & bestlevelPs, double betaM, double betaP, string ofilename);

// function to fit Bayesian density forests under shifts, using Bayesian model averaging
void density_RF(parameters & p, parameters & para, vector<OneDensity_data>& ensemble, int ensemble_sizeP, int n_shiftM, double shift_sizeM, int n_shiftP, double shift_sizeP, double betaM, double betaP, string ofilename, int Cdepth=30);


#endif