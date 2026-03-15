/*
 * File:   SISfunctions.h
 * Author: Mouse
 *
 * Created on April 5, 2012, 1:28 PM
 */

#ifndef SISFUNCTIONS_H
#define	SISFUNCTIONS_H

//#include "hellingerdist.h"
#include "KLdist.h"
#include "data_store.h"

void Create(parameters & p, vector<OnePartition_data>& OrigPs, bool UQ_flag =false);


//void SIS( parameters & p, bool discrete, vector<OnePartition_data>& P_bestlevel, int & bestmodelBP, double beta=1.0);  //commented by thchiu
void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, double beta = 1.0, int Cdepth=30, int max_depth =100, string method="gibbs", double tree_temperature=1.0, bool UQ_flag = false);   //added by thchiu

void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, string ofilename, vector<double> &mmax,vector<double> &mmin,  double beta=1.0, int Cdepth=30, int max_depth=100, string method="gibbs", double tree_temperature=1.0, bool UQ_flag =false);   //added by thchiu

void SISforlinearity( parameters & p, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp);


double Output_distance_2 ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect);

void Output_distance ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect);

void Output_distance_mixnormal ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp,  double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect,vector<OnePartition_data> SepPs );


// The function for calculating weight associated with each forests
vector<double> Gibbs_FOF_weight(vector<vector<double> >& testdata, vector<OneDensity_data> ensemble,  double betaM = 1.0, double betaP= 1.0, double splitP = 0.5, double tree_temperature=1.0, double forest_temperature=1.0);

void GDF_UQ_sampling(vector<OneForest_data> GDF_posterior, vector<OneForest_data> & UQ_samples, int n_sample, double pseudo_count);

void GDF_UQ_path_evaluation(vector<OneForest_data> & UQ_samples, vector<vector<double>> test_points, string out_dir);

#endif	/* SISFUNCTIONS_H */

