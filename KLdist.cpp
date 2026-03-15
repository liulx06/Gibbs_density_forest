#include "KLdist.h"


vector<double> getvariousKLdist( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x),vector<OnePartition_data>& modePs, vector<double>& modewts) {
//    vector<double> HDs;
//    HDs.push_back(KL_fromf(density, p.testdata, modePs[0], p.samplesize, 0.0));
//    HDs.push_back(KL_mix(density, p.testdata, modePs, modewts, p.samplesize, 0.0));
//    HDs.push_back(KL_fromf(density, p.testdata, modePs[0], p.samplesize, smoothneighbordist));
//    HDs.push_back(KL_mix(density, p.testdata, modePs, modewts, p.samplesize, smoothneighbordist));
//    HDs.push_back(KL_fromf(density, p.testdata, modePs[0], p.samplesize, 2 * smoothneighbordist));
//    HDs.push_back(KL_mix(density, p.testdata, modePs, modewts, p.samplesize, 2 * smoothneighbordist));
    return all_KL(density, p.testdata, modePs, modewts, p.samplesize,  smoothneighbordist);
}

vector<double> getvariousKLdist_mixnormal( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp, vector<OnePartition_data>& modePs, vector<double>& modewts,vector<OnePartition_data>& SepPs) {
    vector<double> HDs(6,0);
    HDs[0]=KL_fromf_mixnormal(density,mixp, p.testdata, modePs[0],SepPs, 0.0);
    HDs[1]=KL_fromf_mixnormal_mix(density,mixp, p.testdata, modePs, modewts, p.samplesize,0.0,SepPs);
    HDs[2]=KL_fromf_mixnormal(density,mixp, p.testdata, modePs[0], SepPs, smoothneighbordist);
    HDs[3]=KL_fromf_mixnormal_mix(density,mixp, p.testdata, modePs, modewts, p.samplesize, smoothneighbordist,SepPs);
    HDs[4]=KL_fromf_mixnormal(density,mixp, p.testdata, modePs[0],  SepPs, 2 * smoothneighbordist);
    HDs[5]=KL_fromf_mixnormal_mix(density,mixp, p.testdata, modePs, modewts, p.samplesize, 2 * smoothneighbordist,SepPs);
    return HDs;
}


double KL_fromf_mixnormal_mix(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double>& wts, int samplesize, double lattic, vector<OnePartition_data>& SepPs ) {
    int dim = (int)testdata[0].size();
    double result = 0;
    int k=0;
    vector<vector<double> > Ftestdata;
    if ((int) SepPs.size() > 0)  Ftransform(testdata, Ftestdata, SepPs);
    double d1, d2;
    for (int i = 0; i < (int) testdata.size(); i++) {
        d1 = density(testdata[i], mixp);
        if ((int) SepPs.size() > 0) {
            d2 = densities_from_mixpartition(Ftestdata[i], Ps, wts, samplesize, lattic);
            for (int d = 0; d < dim; d++) {
                vector<double> onepiece(1, testdata[i][d]);
                d2 *= one_density(SepPs[d], onepiece, -1);
            }
        } else d2 = densities_from_mixpartition(testdata[i], Ps, wts, samplesize, lattic);
        if(d2>0 && d1>0 && !isnan(log(d2/d1)))  {
            result+= log(d2/d1);
            k++;
        }
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result =-1.0*result;
    return result;
}


double KL_fromf_mixnormal(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, OnePartition_data &Ps, vector<OnePartition_data>& SepPs,double lattic) {
    double result = 0;
    double d1, d2;
    int dim = (int)testdata[0].size();
    vector<vector<double> > Ftestdata;
     if ((int) SepPs.size() > 0)  Ftransform(testdata, Ftestdata, SepPs);
    int k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {
        d1 = density(testdata[i],mixp);
        if((int)SepPs.size()>0){
            d2 = one_density(Ps, Ftestdata[i], lattic);
            for(int d = 0; d<dim; d++){
                vector<double> onepiece(1,testdata[i][d]);
                d2 *= one_density(SepPs[d], onepiece, -1);
            }
        }
        else d2= one_density(Ps, testdata[i], lattic);
        if(d2>0 && d1>0 && !isnan(log(d2/d1)))  {
            result+= log(d2/d1);
            k++;
        }
    }
//    cout<<"k="<<k<<endl;
    result /=(double)k;
    result =-1.0*result;
    return result;
}

double KL_fromf_AR1(double (*density)(const vector<double> &x, const vector<double> &mu, double rho, double sigma), const vector<double> &mu, double rho, double sigma, vector<vector<double> >& testdata, OnePartition_data &Ps, vector<OnePartition_data>& SepPs) {
    double result = 0;
    double d1, d2;
    int dim = (int)testdata[0].size();
    vector<vector<double> > Ftestdata;
	if ((int) SepPs.size() > 0)  Ftransform(testdata, Ftestdata, SepPs);
	
    int k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {
        double d1 = density(testdata[i], mu, rho, sigma);
        double d2= one_density(Ps, Ftestdata[i], 0.0);
		for(int d=0; d < dim; d++) {
			vector<double> onepiece(1, testdata[i][d]);
			d2 *= one_density(SepPs[d], onepiece, 0.0);
		}
        if(d2>0 && d1>0 && !isnan(log(d2/d1)))  {
            result+= log(d2/d1);
            k++;
        } else if(d2 ==0 && d1 >0){
        	k++;
        } else{
        	cout<<"d2: "<<d2<<" "<<"d1: "<<d1<<" "<<"log(d2/d1): "<<log(d2/d1)<<endl;
        }
		
		if(i < 10) {
			cout<<testdata[i][0]<<" "<<testdata[i][1]<<" "<<d1<<" "<<d2<<endl;
		}
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;
}

double KL_fromf(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, OnePartition_data P, vector<OnePartition_data> marginalP, double lattic1, double lattic2) {
	int testsize;
	testsize = testdata.size();
	vector<vector<double> > trans_test;
	trans_test.resize(testsize);
	int dim;
	dim = testdata[0].size();
	for(int i=0; i < testsize; i++) {
		trans_test[i].resize(dim);
	}
	
	for(int d=0; d< dim; d++) {
		// for KL dist calculation, linxi
		vector<double> onevartest(testsize, 0);
		for(int i=0; i<testsize; i++) {
			onevartest[i] = testdata[i][d];
		}
		// get marginal quantiles for test points, linxi
		vector<double> newtestx = Ftransform(onevartest, marginalP[d]);
		for(int i=0; i< testsize; i++) {
			trans_test[i][d] = newtestx[i];
		}
		
		// checkings, linxi
		cerr<<"orginal test data is: ";
        for(int i=0; i<10; i++) cerr<<onevartest[i]<<" ";
        cerr<<endl;

        cerr<<"transformed test data is: ";
		for(int i=0; i<10; i++) cerr<<newtestx[i]<<" "<<trans_test[i][d]<<" ";
        cerr<<endl;

        vector<double> recoveronevartest = inv_Ftransform(newtestx, marginalP[d]);

		cerr<<"recoverd test data is: ";
        for(int i=0; i<10; i++) cerr<<recoveronevartest[i]<<" ";
        cerr<<endl;
	}
	
    double result = 0;
    double k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {

        double d1 = density(testdata[i]);
		double d2= one_density(P, trans_test[i], lattic1);	
		for(int d=0; d < dim; d++) {
			vector<double> onepiece(1, testdata[i][d]);
			d2 *= one_density(marginalP[d], onepiece, lattic2);
		}
        if(d2>0 && d1>0 && !isnan(log(d2/d1)))  {
            result+= log(d2/d1);
            k++;
        } else if(d2 ==0 && d1 >0){
        	k++;
        } else{
        	cout<<"d2: "<<d2<<" "<<"d1: "<<d1<<" "<<"log(d2/d1): "<<log(d2/d1)<<endl;
        }
		
		if(i < 10) {
			cout<<testdata[i][0]<<" "<<testdata[i][1]<<" "<<d1<<" "<<d2<<endl;
		}
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;
}

double KL_fromf_shift(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OneDensity_data> ensemble){
	int testsize;
	testsize = testdata.size();
	int dim;
	dim = testdata[0].size();
	
	
	vector<double> true_density(testsize,0);
	for(int i=0; i< testsize; i++) {
		if(testdata[i][0] <0.00001 || testdata[i][0]>0.99999) {
			true_density[i]=0;
		}else {
			true_density[i] =density(testdata[i]);			
		}		
	}
	
	vector<double> est_density(testsize, 0);
	for(int l=0; l<(int)ensemble.size(); l++) {		
		vector<vector<double> > trans_test;
		trans_test.clear();
		trans_test.resize(testsize);
		for(int i=0; i < testsize; i++) {
			trans_test[i].resize(dim);
		}
	
		vector<double> onetreeD(testsize, 1.0);
		for(int d=0; d< dim; d++) {
			//vector<double> onevartest(testsize,0);
			vector<double> shifted_data(testsize, 0);
			cerr<<"shift along dimension "<<d<<" is: "<<ensemble[l].shift[d]<<endl;
			
			/*
		    OnePartition_data &op = ensemble[l].marginalP[d];
			
			for (int i = 0; i < (int) op.sregs.size(); i++) {
		        pair<double, double> range;
				range = convert_ranges(op.sregs[i].reg_code[0]);
		        double ub,lb;
		        ub = Ftransform(range.second, ensemble[l].marginalP[d]);
		        // ub = ub * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
				//ub = ub * 2.0 -onefit.shift[d] - 0.00001;
		        lb = Ftransform(range.first, ensemble[l].marginalP[d]);
		        // lb = lb * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
				// lb = lb *2.0 - onefit.shift[d] - 0.00001;
		        //         outfile<<range.first<<" "<<range.second<<" ";
		        cout <<lb << " " << ub <<" "<<endl;
		        cout << range.first << " " << range.second << " "<<endl;
				// calculate the density on a subregion, linxi
		        //cout<< (double)op.sregs[i].num/(double)samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)+(double)p.dim*log(2.0)) <<" "<< op.sregs[i].num<<endl;
		        //cout << (double)op.sregs[i].num/(double)samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)+(double)p.dim*log(2.0)) <<" "<< op.sregs[i].num<<endl;
			}
			*/
			// vector<double> shifted_data(testsize, 0);
			vector<double> newtestx(testsize,0);
			for(int i=0; i< testsize; i++) {
				shifted_data[i] =0.5*(testdata[i][d] + ensemble[l].shift[d]);
				// newtestx[i] = Ftransform(shifted_data[i], ensemble[l].marginalP[d]);
				// get marginal quantiles for test points based on lth estimate, linxi
				// trans_test[i][d] = 0.5*(newtestx[i]+ensemble[l].shift[d]);
				
				//if(i-i/100==0) {
				//	cout<<i<<" ";
				//}
			}
			
			newtestx = Ftransform(shifted_data, ensemble[l].marginalP[d]);
			
			for(int i=0; i< testsize; i++) {
				// get marginal quantiles for test points based on lth estimate, linxi
				trans_test[i][d] = 0.5*(newtestx[i]+ensemble[l].joint_shift[d]);
			}
			
			// checkings, linxi
			/*
			cerr<<"orginal test data is: ";
	        for(int i=0; i<10; i++) cerr<<shifted_data[i]<<" "<<testdata[i][d]<<" ";
	        cerr<<endl;

	        cerr<<"transformed test data is: ";
			for(int i=0; i<10; i++) cerr<<newtestx[i]<<" "<<trans_test[i][d]<<" ";
	        cerr<<endl;

	        vector<double> recoveronevartest = inv_Ftransform(newtestx, ensemble[l].marginalP[d]);

			cerr<<"recoverd test data is: ";
	        for(int i=0; i<10; i++) cerr<<recoveronevartest[i]<<" "<<2*recoveronevartest[i]-ensemble[l].shift[d]<<" ";
	        cerr<<endl;
			*/
			
			// calculate estimated density based on one tree (product of marginal densities)
			for(int i=0; i < testsize; i++) {
				vector<double> onepiece(1, shifted_data[i]);
				onetreeD[i] *= one_density(ensemble[l].marginalP[d], onepiece, 0.0);
			}
			
		}
		
		for(int i=0; i< testsize; i++) {
			if(dim > 1) {
				// calculate estimated density based on one tree (joint density)
				onetreeD[i]*= one_density(ensemble[l].jointP, trans_test[i],0.0);
				onetreeD[i]*= exp(-(double)(2*dim)*log(2.0));
			
				// update estimated density based on the ensemble
				est_density[i] += onetreeD[i];
			}
			else{
				onetreeD[i]*= exp(-(double)(dim)*log(2.0));
				est_density[i] += onetreeD[i];
			}
		}		
	
	}
	
	for(int i=0; i < testsize; i++) {
		if(testdata[i][0] <0.00001 || testdata[i][0]>0.99999) {
			est_density[i]=0;
		}else {
			est_density[i] *= 1/(double)ensemble.size();
		}		
	}
	
    double k=0;
	double result =0;
    for (int i = 0; i < testsize; i++) {
		
        if(true_density[i]>0 && est_density[i]>0 && !isnan(log(est_density[i]/true_density[i])))  {
            result+= log(est_density[i]/true_density[i]);
            k++;
        } else if(est_density[i]==0 &&  true_density[i] >0){
        	k++;
        } else{
        	cout<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<" "<<"log(est/true): "<<log(est_density[i]/true_density[i])<<endl;
        }
		
		if(i < 10) {
			cout<<testdata[i][0]<<" "<<testdata[i][1]<<" "<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<endl;
		}
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;

}

double KL_fromf_shift_AR1(double (*density)(const vector<double> &x, const vector<double> &mu, double rho, double sigma), const vector<double> &mu, double rho, double sigma, vector<vector<double> >& testdata, vector<OneDensity_data> ensemble){
	int testsize;
	testsize = testdata.size();
	int dim;
	dim = testdata[0].size();
	
	
	vector<double> true_density(testsize,0);
	for(int i=0; i< testsize; i++) {
		true_density[i] =density(testdata[i], mu, rho, sigma);
	}
	
	vector<double> est_density(testsize, 0);
	for(int l=0; l<(int)ensemble.size(); l++) {		
		vector<vector<double> > trans_test;
		trans_test.clear();
		trans_test.resize(testsize);
		for(int i=0; i < testsize; i++) {
			trans_test[i].resize(dim);
		}
	
		vector<double> onetreeD(testsize, 1.0);
		for(int d=0; d< dim; d++) {
			//vector<double> onevartest(testsize,0);
			vector<double> shifted_data(testsize, 0);
			cerr<<"shift along dimension "<<d<<" is: "<<ensemble[l].shift[d]<<endl;
			
			vector<double> newtestx(testsize,0);
			for(int i=0; i< testsize; i++) {
				shifted_data[i] =0.5*(testdata[i][d] + ensemble[l].shift[d]);
			}
			
			newtestx = Ftransform(shifted_data, ensemble[l].marginalP[d]);
			
			for(int i=0; i< testsize; i++) {
				// get marginal quantiles for test points based on lth estimate, linxi
				trans_test[i][d] = 0.5*(newtestx[i]+ensemble[l].joint_shift[d]);
			}
			
			for(int i=0; i < testsize; i++) {
				vector<double> onepiece(1, shifted_data[i]);
				onetreeD[i] *= one_density(ensemble[l].marginalP[d], onepiece, 0.0);
			}
			
		}
		
		for(int i=0; i< testsize; i++) {
			if(dim > 1) {
				// calculate estimated density based on one tree (joint density)
				onetreeD[i]*= one_density(ensemble[l].jointP, trans_test[i],0.0);
				onetreeD[i]*= exp(-(double)(2*dim)*log(2.0));
			
				// update estimated density based on the ensemble
				est_density[i] += onetreeD[i];
			}
			else{
				onetreeD[i]*= exp(-(double)dim *log(2.0));
				est_density[i] += onetreeD[i];
			}
		}		
	
	}
	
	for(int i=0; i < testsize; i++) {
		est_density[i] *= 1/(double)ensemble.size();
	}
	
    double k=0;
	double result =0;
    for (int i = 0; i < testsize; i++) {
		
        if(true_density[i]>0 && est_density[i]>0 && !isnan(log(est_density[i]/true_density[i])))  {
            result+= log(est_density[i]/true_density[i]);
            k++;
        } else if(est_density[i]==0 &&  true_density[i] >0){
        	k++;
        } else{
        	cout<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<" "<<"log(est/true): "<<log(est_density[i]/true_density[i])<<endl;
        }
		
		if(i < 10) {
			cout<<testdata[i][0]<<" "<<testdata[i][1]<<" "<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<endl;
		}
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;

}

double KL_fromf_shift_mix2Dnorm(double (*density)(const vector<double> &x, mixnormalparas & mixp), mixnormalparas & mixp, vector<vector<double> >& testdata, vector<OneDensity_data> ensemble){
	int testsize;
	testsize = testdata.size();
	int dim;
	dim = testdata[0].size();
	
	
	vector<double> true_density(testsize,0);
	for(int i=0; i< testsize; i++) {
		true_density[i] =density(testdata[i], mixp);
	}
	
	vector<double> est_density(testsize, 0);
	for(int l=0; l<(int)ensemble.size(); l++) {		
		vector<vector<double> > trans_test;
		trans_test.clear();
		trans_test.resize(testsize);
		for(int i=0; i < testsize; i++) {
			trans_test[i].resize(dim);
		}
	
		vector<double> onetreeD(testsize, 1.0);
		for(int d=0; d< dim; d++) {
			//vector<double> onevartest(testsize,0);
			vector<double> shifted_data(testsize, 0);
			cerr<<"shift along dimension "<<d<<" is: "<<ensemble[l].shift[d]<<endl;
			
			vector<double> newtestx(testsize,0);
			for(int i=0; i< testsize; i++) {
				shifted_data[i] =0.5*(testdata[i][d] + ensemble[l].shift[d]);
			}
			
			newtestx = Ftransform(shifted_data, ensemble[l].marginalP[d]);
			
			for(int i=0; i< testsize; i++) {
				// get marginal quantiles for test points based on lth estimate, linxi
				trans_test[i][d] = 0.5*(newtestx[i]+ensemble[l].joint_shift[d]);
			}
			
			for(int i=0; i < testsize; i++) {
				vector<double> onepiece(1, shifted_data[i]);
				onetreeD[i] *= one_density(ensemble[l].marginalP[d], onepiece, 0.0);
			}
			
		}
		
		for(int i=0; i< testsize; i++) {
			if(dim > 1) {
				// calculate estimated density based on one tree (joint density)
				onetreeD[i]*= one_density(ensemble[l].jointP, trans_test[i],0.0);
				onetreeD[i]*= exp(-(double)(2*dim)*log(2.0));
			
				// update estimated density based on the ensemble
				est_density[i] += onetreeD[i];
			}
			else{
				onetreeD[i]*= exp(-(double)dim *log(2.0));
				est_density[i] += onetreeD[i];
			}
		}		
	
	}
	
	for(int i=0; i < testsize; i++) {
		est_density[i] *= 1/(double)ensemble.size();
	}
	
    double k=0;
	double result =0;
    for (int i = 0; i < testsize; i++) {
		
        if(true_density[i]>0 && est_density[i]>0 && !isnan(log(est_density[i]/true_density[i])))  {
            result+= log(est_density[i]/true_density[i]);
            k++;
        } else if(est_density[i]==0 &&  true_density[i] >0){
        	k++;
        } else{
        	cout<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<" "<<"log(est/true): "<<log(est_density[i]/true_density[i])<<endl;
        }
		
		if(i < 10) {
			cout<<testdata[i][0]<<" "<<testdata[i][1]<<" "<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<endl;
		}
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;

}

double KL_fromf_fof(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OneForest_data> forest_of_Fs){
	int testsize;
	testsize = testdata.size();
	int dim;
	dim = testdata[0].size();
	
	int n_forests =(int)forest_of_Fs.size();
	vector<double> lweights(n_forests, 0), weights(n_forests, 0);
	double w_lsum =0.0;
	for(int k=0; k<n_forests; k++) {
		lweights[k] = forest_of_Fs[k].weight;
		cout<<"weight associated with forest "<<k<<" is:"<<forest_of_Fs[k].weight<<endl;
		w_lsum += lweights[k];
	}
	w_lsum = w_lsum/(double)n_forests;
	cout<<"w_lsum is: "<<w_lsum<<endl;
	if(isnan(w_lsum)) {
		for(int k=0; k< n_forests; k++) {
			lweights[k] =0;
			cout<<"asjusted weight associated with forest "<<k<<" is:"<<lweights[k]<<endl;
		}
	}else{
		for(int k=0; k<n_forests; k++) {
			lweights[k] = lweights[k]- w_lsum;
			cout<<"asjusted weight associated with forest "<<k<<" is:"<<lweights[k]<<endl;
		}	
	}	
	w_lsum = logsum(lweights);
	cout<<"w_lsum is: "<<w_lsum<<endl;
	double w_sum =0.0;
	for(int k=0; k<n_forests; k++) {
		lweights[k] -= w_lsum;
		weights[k] = exp(lweights[k]);
		w_sum += weights[k];
	}
	cout<<"summation of normalizeds weights is: "<<w_sum<<endl;
	
	vector<double> true_density(testsize,0);
	for(int i=0; i< testsize; i++) {
		if(testdata[i][0] <0.00001 || testdata[i][0]>0.99999) {
			true_density[i]=0;
		}else {
			true_density[i] =density(testdata[i]);			
		}		
	}
	
	vector<double> est_density(testsize, 0);
	for(int k=0; k<n_forests; k++) {
		vector<double> oneforestD(testsize,0);
		for(int l=0; l<(int)forest_of_Fs[k].forestP.size(); l++) {		
			vector<vector<double> > trans_test;
			trans_test.clear();
			trans_test.resize(testsize);
			for(int i=0; i < testsize; i++) {
				trans_test[i].resize(dim);
			}
	
			vector<double> onetreeD(testsize, 1.0);
			for(int d=0; d< dim; d++) {
				//vector<double> onevartest(testsize,0);
				vector<double> shifted_data(testsize, 0);
				cerr<<"shift along dimension "<<d<<" is: "<<forest_of_Fs[k].forestP[l].shift[d]<<endl;
			
				/*
			    OnePartition_data &op = ensemble[l].marginalP[d];
			
				for (int i = 0; i < (int) op.sregs.size(); i++) {
			        pair<double, double> range;
					range = convert_ranges(op.sregs[i].reg_code[0]);
			        double ub,lb;
			        ub = Ftransform(range.second, ensemble[l].marginalP[d]);
			        // ub = ub * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
					//ub = ub * 2.0 -onefit.shift[d] - 0.00001;
			        lb = Ftransform(range.first, ensemble[l].marginalP[d]);
			        // lb = lb * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
					// lb = lb *2.0 - onefit.shift[d] - 0.00001;
			        //         outfile<<range.first<<" "<<range.second<<" ";
			        cout <<lb << " " << ub <<" "<<endl;
			        cout << range.first << " " << range.second << " "<<endl;
					// calculate the density on a subregion, linxi
			        //cout<< (double)op.sregs[i].num/(double)samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)+(double)p.dim*log(2.0)) <<" "<< op.sregs[i].num<<endl;
			        //cout << (double)op.sregs[i].num/(double)samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)+(double)p.dim*log(2.0)) <<" "<< op.sregs[i].num<<endl;
				}
				*/
				// vector<double> shifted_data(testsize, 0);
				vector<double> newtestx(testsize,0);
				for(int i=0; i< testsize; i++) {
					shifted_data[i] =0.5*(testdata[i][d] + forest_of_Fs[k].forestP[l].shift[d]);
					// newtestx[i] = Ftransform(shifted_data[i], ensemble[l].marginalP[d]);
					// get marginal quantiles for test points based on lth estimate, linxi
					// trans_test[i][d] = 0.5*(newtestx[i]+ensemble[l].shift[d]);
				
					//if(i-i/100==0) {
					//	cout<<i<<" ";
					//}
				}
			
				newtestx = Ftransform(shifted_data, forest_of_Fs[k].forestP[l].marginalP[d]);
			
				for(int i=0; i< testsize; i++) {
					// get marginal quantiles for test points based on lth estimate, linxi
					trans_test[i][d] = 0.5*(newtestx[i]+forest_of_Fs[k].forestP[l].joint_shift[d]);
				}
			
				// checkings, linxi
				/*
				cerr<<"orginal test data is: ";
		        for(int i=0; i<10; i++) cerr<<shifted_data[i]<<" "<<testdata[i][d]<<" ";
		        cerr<<endl;

		        cerr<<"transformed test data is: ";
				for(int i=0; i<10; i++) cerr<<newtestx[i]<<" "<<trans_test[i][d]<<" ";
		        cerr<<endl;

		        vector<double> recoveronevartest = inv_Ftransform(newtestx, ensemble[l].marginalP[d]);

				cerr<<"recoverd test data is: ";
		        for(int i=0; i<10; i++) cerr<<recoveronevartest[i]<<" "<<2*recoveronevartest[i]-ensemble[l].shift[d]<<" ";
		        cerr<<endl;
				*/
			
				// calculate estimated density based on one tree (product of marginal densities)
				for(int i=0; i < testsize; i++) {
					vector<double> onepiece(1, shifted_data[i]);
					onetreeD[i] *= one_density(forest_of_Fs[k].forestP[l].marginalP[d], onepiece, 0.0);
				}
			
			}
		
			for(int i=0; i< testsize; i++) {
				if(dim > 1) {
					// calculate estimated density based on one tree (joint density)
					onetreeD[i]*= one_density(forest_of_Fs[k].forestP[l].jointP, trans_test[i],0.0);
					onetreeD[i]*= exp(-(double)(2*dim)*log(2.0));
			
					// update estimated density based on the ensemble
					oneforestD[i] += onetreeD[i];
				}
				else{
					onetreeD[i]*= exp(-(double)(dim)*log(2.0));
					oneforestD[i] += onetreeD[i];
				}
			}		
	
		}	
		for(int i=0; i<testsize; i++) {
			oneforestD[i] *= 1/(double)forest_of_Fs[k].forestP.size();
			est_density[i] += oneforestD[i] *weights[k];
		}
	}
	
	for(int i=0; i < testsize; i++) {
		if(testdata[i][0] <0.00001 || testdata[i][0]>0.99999) {
			est_density[i]=0;
		}
	}
	
    int count=0;
	double result =0;
    for (int i = 0; i < testsize; i++) {
		
        if(true_density[i]>0 && est_density[i]>0 && !isnan(log(est_density[i]/true_density[i])))  {
            result+= log(est_density[i]/true_density[i]);
            count++;
        } else if(est_density[i]==0 &&  true_density[i] >0){
        	count++;
        } else{
        	cout<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<" "<<"log(est/true): "<<log(est_density[i]/true_density[i])<<endl;
        }
		
		if(i < 10) {
			cout<<testdata[i][0]<<" "<<testdata[i][1]<<" "<<"estimated density: "<<est_density[i]<<" "<<"true density: "<<true_density[i]<<endl;
		}
    }
    cout<<"count="<<count<<endl;
    result /=(double)count;
    result = -1.0*result;
    return result;

}

double KL_mix(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic) {

    double result = 0;
     double k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {

        double d1 = density(testdata[i]);
        double d2 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, lattic);
        if(d2>0 && d1>0&& !isnan(log(d2/d1))) {
             k++;
             result+= log(d2/d1);
        }
    }

    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;
}

// the best partition should be Ps[0]
// potentially can be faster if combine lattice=0 and >0
vector<double> all_KL(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic) {

    vector<double> dists(6,0);
    double k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {

        double d1 = density(testdata[i]);
        double d2 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, lattic);
        double d3 = one_density(Ps[0], testdata[i], lattic);
        double d4 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, 0.0);
        double d5 = one_density(Ps[0], testdata[i], 0.0);

        double d6 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, lattic*2);
        double d7 = one_density(Ps[0], testdata[i], lattic*2);
        if(d2>0 && d1>0 && d3>0 && d4>0 && d5>0 && d6>0&& d7>0&&!isnan(log(d2/d1))) {
             k++;
             dists[0]+= log(d5/d1);
             dists[1]+= log(d4/d1);
             dists[2]+= log(d3/d1);
             dists[3]+= log(d2/d1);
             dists[4]+= log(d7/d1);
             dists[5]+= log(d6/d1);
        }
    }

    cout<<"k="<<k<<endl;
    for(int i=0; i<6; i++) {
        dists[i] /=(double)k;
        dists[i]= -dists[i];
    }
    return dists;
}


double KL_fromf_mixnormal_compare(double (*density)(const vector<double> &x, mixnormalparas &mixp),mixnormalparas &mixp, vector<vector<double> >& testdata, vector<double>  matlabresult) {

    double result = 0;
    int k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {
        double d1 = density(testdata[i], mixp);
        double d2 = matlabresult[i];
        if(d2>0 && d1>0&& !isnan(log(d2/d1))) {
             k++;
             result+= log(d2/d1);
        }
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result =-1.0*result;
    return result;
}

double KL_fromf_compare(double (*density)(const vector<double> &x), vector<vector<double> >& testdata, vector<double>  matlabresult) {

    double result = 0;
    int k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {
        double d1 = density(testdata[i]);
        double d2 = matlabresult[i];
        if(d2>0 && d1>0&& !isnan(log(d2/d1))) {
             k++;
             result+= log(d2/d1);
        }
    }
    result /=(double)k;
    result =-1.0*result;
    return result;
}

// Calculate the KL for different settings in simulation studies, Linxi
//double simulation_KL_fromf(vector<vector<double> >& testdata, OnePartition_data P, vector<OnePartition_data> marginalP);

double simulation_KL_fromf_shift(vector<vector<double> >& testdata, vector<OneDensity_data> ensemble, string data_model) {
	double KLdist_shift =0; 
	if(data_model.find("normAR1") != string::npos) {		
		size_t coef_pos = data_model.find(".");
		string modelCoef = data_model.substr(coef_pos);
		double rho = strTo<double>(modelCoef);
		cout<<"autocorrelation coefficient is: "<<modelCoef<<" "<<strTo<double>(modelCoef)<<endl;
		
		int dim = (int)testdata[0].size();
		vector<double> normmu(dim, 0.5);
		KLdist_shift = KL_fromf_shift_AR1(sample_normal_AR1_density, normmu, rho, 0.1, testdata, ensemble);	
		return KLdist_shift; 
	} else if(data_model =="mix2Dnorm") {
		mixnormalparas mixPara;
		mixPara.mu1 = 0.25;
		mixPara.mu2 = 0.25;
		mixPara.mu3 = 0.75;
		mixPara.mu4 = 0.75;
		mixPara.sd1 = 0.05;
		mixPara.sd2 = 0.05;
		mixPara.ratio = 0.4;
	
		KLdist_shift = KL_fromf_shift_mix2Dnorm(sample_mix2normal_density, mixPara, testdata, ensemble);
		return KLdist_shift; 
	} else if(data_model =="1Dbeta") {
		KLdist_shift = KL_fromf_shift(sample_beta15_density, testdata, ensemble);
		return KLdist_shift; 
	} else if(data_model == "1Dcos") {
		KLdist_shift = KL_fromf_shift(sample_cos_density, testdata, ensemble);
		return KLdist_shift; 
	}else if (data_model == "1Dcircle") {
		KLdist_shift = KL_fromf_shift(sample_semi_circle_density, testdata, ensemble);
		return KLdist_shift; 
	}else if(data_model.find("mix5Dnorm") != string::npos) {
		KLdist_shift = KL_fromf_shift(sample_mix5Dnormal_density, testdata, ensemble);
		return KLdist_shift; 
	} else {
		cout<<"The model for data generating is not found"<<endl;
		return -1.0;
	}
}

vector<double> simulation_KL_fromf_fof(vector<vector<double> >& testdata, vector<OneForest_data> forest_of_Fs, string data_model, int bestFind) {	
	double KLdist_Gibbs=0.0; 
	// calculating KL divergence, linxi
	double KLdist_fof=0.0;
	if(data_model.find("normAR1") != string::npos) {		
		size_t coef_pos = data_model.find(".");
		string modelCoef = data_model.substr(coef_pos);
		double rho = strTo<double>(modelCoef);
		cout<<"autocorrelation coefficient is: "<<modelCoef<<" "<<strTo<double>(modelCoef)<<endl;
	
		// vector<double> normmu(p.dim, 0.5);
		//KLdist_shift = KL_fromf_shift_AR1(sample_normal_AR1_density, normmu, rho, 0.1, p.testdata, ensemble);	
	} else if(data_model =="mix2Dnorm") {
		mixnormalparas mixPara;
		mixPara.mu1 = 0.25;
		mixPara.mu2 = 0.25;
		mixPara.mu3 = 0.75;
		mixPara.mu4 = 0.75;
		mixPara.sd1 = 0.05;
		mixPara.sd2 = 0.05;
		mixPara.ratio = 0.4;
	
		//KLdist_shift = KL_fromf_shift_mix2Dnorm(sample_mix2normal_density, mixPara, p.testdata, ensemble);
	} else if(data_model =="1Dbeta") {
		KLdist_fof = KL_fromf_fof(sample_beta15_density, testdata, forest_of_Fs);
		KLdist_Gibbs = KL_fromf_shift(sample_beta15_density, testdata, forest_of_Fs[bestFind].forestP);
	} else if(data_model == "1Dcos") {
		KLdist_fof = KL_fromf_fof(sample_cos_density, testdata, forest_of_Fs);
		KLdist_Gibbs = KL_fromf_shift(sample_cos_density, testdata, forest_of_Fs[bestFind].forestP);
	}else if (data_model == "1Dcircle") {
		KLdist_fof = KL_fromf_fof(sample_semi_circle_density, testdata, forest_of_Fs);
		KLdist_Gibbs = KL_fromf_shift(sample_semi_circle_density, testdata, forest_of_Fs[bestFind].forestP);
	}else if(data_model.find("mix5Dnorm") != string::npos) {
		KLdist_fof = KL_fromf_fof(sample_mix5Dnormal_density, testdata, forest_of_Fs);
		KLdist_Gibbs = KL_fromf_shift(sample_mix5Dnormal_density, testdata, forest_of_Fs[bestFind].forestP);
	} else {
		cout<<"The model for data generating is not found"<<endl;
		KLdist_fof = -1.0;
		KLdist_Gibbs = -1.0;
	}
	
	vector<double> KLdist_GDF = {0.0, 0.0};
	KLdist_GDF[0] = KLdist_Gibbs; // KL for posterior mode
	KLdist_GDF[1] = KLdist_fof; // KL for posterior mean
	return KLdist_GDF;
}

double simulation_KL_fromf(vector<vector<double> >& testdata, OnePartition_data bestlevelPs, vector<OnePartition_data> marginalP, string data_model) {
	double KLdist=0.0;
	if(data_model.find("normAR1") != string::npos) {		
		size_t coef_pos = data_model.find(".");
		string modelCoef = data_model.substr(coef_pos);
		double rho = strTo<double>(modelCoef);
		cout<<"autocorrelation coefficient is: "<<modelCoef<<" "<<strTo<double>(modelCoef)<<endl;
		
		int dim = (int)testdata[0].size();
		vector<double> normmu(dim, 0.5);
		KLdist = KL_fromf_AR1(sample_normal_AR1_density, normmu, rho, 0.1, testdata, bestlevelPs, marginalP);
		//KLdist = KL_fromf_shift_AR1(sample_normal_AR1_density, normmu, rho, 0.1, p.testdata, ensemble);	
	} else if(data_model =="mix2Dnorm") {
		mixnormalparas mixPara;
		mixPara.mu1 = 0.25;
		mixPara.mu2 = 0.25;
		mixPara.mu3 = 0.75;
		mixPara.mu4 = 0.75;
		mixPara.sd1 = 0.05;
		mixPara.sd2 = 0.05;
		mixPara.ratio = 0.4;
	
		KLdist = KL_fromf_mixnormal(sample_mix2normal_density, mixPara, testdata, bestlevelPs, marginalP, 0.0);
		//KLdist = KL_fromf_shift_mix2Dnorm(sample_mix2normal_density, mixPara, p.testdata, ensemble);
	} else if(data_model =="1Dbeta") {
		KLdist = KL_fromf(sample_beta15_density, testdata, bestlevelPs, marginalP, 0.0, 0.0);
	} else if(data_model == "1Dcos") {
		KLdist = KL_fromf(sample_cos_density, testdata, bestlevelPs, marginalP, 0.0, 0.0);
	} else if (data_model == "1Dcircle") {
		KLdist = KL_fromf(sample_semi_circle_density, testdata, bestlevelPs, marginalP, 0.0, 0.0);
	}else if(data_model.find("mix5Dnorm") != string::npos) {
		KLdist = KL_fromf(sample_mix5Dnormal_density, testdata, bestlevelPs, marginalP, 0.0, 0.0);
	} else {
		cout<<"The model for data generating is not found"<<endl;
	}
	
	return KLdist;
}
