#ifndef RANDOMSHIFT_H
#define	RANDOMSHIFT_H

#include "sampling.h"
#include "output.h"


void scale_data_for_shift(vector<vector<double> > & data, double s);

void shift(vector<vector<double> > & data, double h);