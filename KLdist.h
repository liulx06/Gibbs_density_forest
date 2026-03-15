/*
 * File:   KLdist.h
 * Author: Mouse
 *
 * Created on May 21, 2012, 8:21 PM
 */

#ifndef KLDIST_H
#define	KLDIST_H

#include "sampling.h"
#include "output.h"
#include "sample_density_generation.h"
//#include "hellingerdist.h"


vector<double> getvariousKLdist( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x),vector<OnePartition_data>& modePs, vector<double>& modewts);

vector<double> getvariousKLdist_mixnormal( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp, vector<OnePartition_data>& modePs, vector<double>& modewts,vector<OnePartition_data>& SepPs);

vector<double> all_KL(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic);

double KL_fromf_mixnormal_mix(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double>& wts, int samplesize, double lattic, vector<OnePartition_data>& SepPs );

double KL_fromf_mixnormal(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, OnePartition_data &Ps, vector<OnePartition_data>& SepPs, double lattic);

double KL_fromf_AR1(double (*density)(const vector<double> &x, const vector<double> &mu, double rho, double sigma), const vector<double> &mu, double rho, double sigma, vector<vector<double> >& testdata, OnePartition_data &Ps, vector<OnePartition_data>& SepPs);
	
double KL_fromf(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, OnePartition_data P, vector<OnePartition_data> marginalP, double lattic1, double lattic2);

double KL_fromf_shift(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OneDensity_data> ensemble);

double KL_fromf_shift_AR1(double (*density)(const vector<double> &x, const vector<double> &mu, double rho, double sigma), const vector<double> &mu, double rho, double sigma, vector<vector<double> >& testdata, vector<OneDensity_data> ensemble);

double KL_fromf_shift_mix2Dnorm(double (*density)(const vector<double> &x, mixnormalparas & mixp), mixnormalparas & mixp, vector<vector<double> >& testdata, vector<OneDensity_data> ensemble);

double KL_fromf_fof(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OneForest_data> forest_of_Fs);	

double KL_mix(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic);

double KL_fromf_mixnormal_compare(double (*density)(const vector<double> &x, mixnormalparas &mixp),mixnormalparas &mixp, vector<vector<double> >& testdata, vector<double>  matlabresult);

double KL_fromf_compare(double (*density)(const vector<double> &x), vector<vector<double> >& testdata, vector<double>  matlabresult);

// Calculate the KL for different settings in simulation studies, Linxi
double simulation_KL_fromf(vector<vector<double> >& testdata, OnePartition_data bestlevelPs, vector<OnePartition_data> marginalP, string data_model);

double simulation_KL_fromf_shift(vector<vector<double> >& testdata, vector<OneDensity_data> ensemble, string data_model);

vector<double> simulation_KL_fromf_fof(vector<vector<double> >& testdata, vector<OneForest_data> forest_of_Fs, string data_model, int bestFind);

#endif	/* KLDIST_H */

