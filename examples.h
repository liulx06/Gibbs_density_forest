/*
 * File:   examples.h
 * Author: Mouse
 *
 * Created on May 21, 2012, 3:14 PM
 */

#ifndef EXAMPLES_H
#define	EXAMPLES_H

#include "sample_density_generation.h"
#include "SISfunctions.h"
#include "buildpathcttable.h"

void FSIS2(parameters & p,bool knowndensity, double (*density)(const vector<double> &x),   bool  pathctcorrect, bool discrete, double smoothneighbordist, vector< map< OnePartition , double, CompairSRegs> >& pathctmap);
void cytometry();
void PK();
void cash();

#endif	/* EXAMPLES_H */

