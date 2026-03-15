
#include "SISfunctions.h"
#include "timer_ed520.h"    //+ by ed520
extern TIMER timer;

void Create(parameters & p, vector<OnePartition_data>& OrigPs, bool UQ_flag){
    usint_mask um;
    um.mask = 0;
    um.x = 0;
    vector<usint_mask> root(p.dim, um);    //(0,0) should be the root of the sample space.
    shrink_reg_data_short sr;
    sr.num = p.data.size();
	sr.valiNum = p.validata.size();
	if(UQ_flag) {
		sr.inferNum = p.inferdata.size();		
	}	
	sr.regionID =0;
    //data_store<vector<double> > originaldata(p.data); //- by ed520
    //sr.regdata = originaldata;    //- by ed520
    sr.reg_code = root;
    OnePartition_data Pt;
    Pt.w = 0;
	Pt.BF =0.0; // This is initialized as the uniform distribution on the unit cube
	Pt.depth =0;
	Pt.prior = 0.0; // -log(2.0)
	Pt.post_p = 0.0; // - log(2.0) This is proportional to the posterior probability of the root, prob to make a split = 0.5 * 2^{- \nu *depth}
    Pt.sregs.push_back(sr);
	Pt.cut_index.clear();
	Pt.left_count.clear();
	Pt.right_count.clear();
    for (int i = 0; i < p.n; i++) {
        Pt.sregs[0].region_start = i * p.samplesize;    //+ by thchiu + by ed520
		Pt.sregs[0].region_start_vali = i*p.valisize;
		if(p.dim ==1 && UQ_flag) {
			Pt.sregs[0].region_start_infer = i*p.infersize;
		}
        Pt.PartitionID=i;   //+ by thchiu
        OrigPs.push_back(Pt);
    }
 }

bool is_milestone(int i, int plevels){
    if(i==(int)((double)plevels/10.0)) return true;
    if(i==(int)((double)plevels/5.0)) return true;
    if(i==(int)((double)plevels/3.0)) return true;
    if(i==(int)((double)plevels/2.0)) return true;
    if(i==(int)((double)plevels*2.0/3.0)) return true;
    if(i==(int)((double)plevels*4.0/5.0)) return true;
    return false;
}

void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, double beta, int Cdepth, int max_depth, string method, double tree_temperature, bool UQ_flag)   //added by thchiu
{
    vector<double> mmax = vector<double>();
    vector<double> mmin = vector<double>();

    SIS(p, discrete, P_bestlevel, "", mmax, mmin, beta, Cdepth, max_depth, method, tree_temperature, UQ_flag);
}
//void SIS( parameters & p, bool discrete, vector<OnePartition_data>& P_bestlevel, int & bestmodelBP, double beta) {  //commented by thchiu
void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, string ofilename, vector<double>& mmax,vector<double> &mmin, double beta, int Cdepth, int max_depth, string method, double tree_temperature, bool UQ_flag) { //added by thchiu
#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerStart_cpu(timer.sTime_C_SIS);
    else
        timer.TimerStart_cpu(timer.sTime_NC_SIS);
#endif
    cerr<<"beta="<<beta<<endl;
    double* time_array = NULL;
    vector<vector<double> >* c_originals = new vector<vector<double> >;   //+ by ed520
    (*c_originals).resize(p.n);  //+ by ed520
    vector<vector<double> >* c_originalsPtr;    //+ by thchiu
    c_originalsPtr = c_originals;  //+ by thchiu
    vector<int>* cutregs = new vector<int>;    //the region to be cut, thchiu + by ed520
    (*cutregs).resize(p.n);  //+ by ed520
    vector<int>* cutregsPtr; //+ by thchiu
    cutregsPtr = cutregs;  //+ by thchiu
    vector<OnePartition_data> *newPsPtr = new vector<OnePartition_data>; // + by ed520
    Create(p, *newPsPtr, UQ_flag);       //+ by ed520
    vector<OnePartition_data>* PsPtr;   //added by thchiu
    PsPtr = newPsPtr;     //added by thchiu + by ed520

    int modelBP=0;
    double currlBPscore = -1.0; // negative one means no cut...


    // added by johnmu: need the case when you don't cut
    P_bestlevel = (*PsPtr)[0];

    double maxlBPscore=currlBPscore;
    int maxlevel=0;
	
	// For random forest, the stopping rule is determined by tree depth
	if(method == "RF") {
		double tempNum = exp((double)(Cdepth) * log(2.0));
		p.levels = (int)tempNum;
		if((double)p.levels +1 - tempNum < tempNum - (double)p.levels) p.levels +=1;
	}

    for (int i = 0; i < p.levels; i++) {
        cout << "\n\n#i=" << i << endl;
        if ((i > 0) && (p.resampling) && (i % p.steps == 0) && (method != "RF")) { //(ess<10) &&  with getmax included.
            resample(PsPtr, p, c_originalsPtr, cutregsPtr, UQ_flag); //added by thchiu
        }
        bool exhaust = sample_trees(*PsPtr, *c_originalsPtr, *cutregsPtr, p.maxpercentage, time_array, discrete, p, beta, Cdepth, method, tree_temperature, UQ_flag);  //+ by thchiu
        
        if (exhaust) break; // in case there are no more regions to split(letter)
		
		/*
		int mdepth =0;
		
		for(int l=0; l< p.n; l++) {
			if((*PsPtr)[l].depth > mdepth) mdepth = (*PsPtr)[l].depth; 
		}
		
		if(mdepth >= max_depth) {
		    vector<double> ec(p.n, 0);
		    for (int l = 0; l < p.n; l++) {
				ec[l] = exp((*PsPtr)[l].w);
		    }
			int index = rand_int(ec);
			P_bestlevel = (*PsPtr)[index];
            cout<<"the max depth "<<max_depth<<"is achieved!\n";
            break;
		}
		*/
        
        
        //currlBPscore = -lgamma_c(p.samplesize + (i + 2) / 2.0)- (i + 2) * lgamma_c(0.5) + lgamma_c(0.5 * (i + 2)) - beta * (i + 2);
        
        if(method != "RF"){
	        currlBPscore = getmaxlBP(*PsPtr, modelBP);    //+ by thchiu
	        cout << currlBPscore << "\n";

	        if (currlBPscore > maxlBPscore) {
	            maxlBPscore = currlBPscore;
	            //P_bestlevel = newPs;    //problem4: data copy, I think only return the best partition is OK, return the whole partitions isn't effciency, thchiu, - by thchiu
	            P_bestlevel = (*PsPtr)[modelBP];//+ by thchiu
	            cout<<"modelBP="<<modelBP<<endl;
	            //bestmodelBP = modelBP;  //problem4:If we return the best partition, we don't need this anymore. commented by thchiu
	            maxlevel=i;
	            cout<<"newmax level="<<i+2<<endl;
	        }
			// print intermediate results, linxi
	        if( is_milestone(i, p.levels) && ofilename!="" ) {
	             string ofilename2 =ofilename+"_level_"+toStr<int>(i)+".txt";
	             ofstream outfile (ofilename2.c_str());
	             print_partition(outfile,P_bestlevel, mmax, mmin);
	        }
	        if(i>maxlevel+10) {
				/*
			    vector<double> ec(p.n, 0);
			    for (int l = 0; l < p.n; l++) {
					ec[l] = exp((*PsPtr)[l].w);
			    }
				int index = rand_int(ec);
				P_bestlevel = (*PsPtr)[index];
	            cout<<"no improvement in 10 levels!\n";
	            break;
				*/
				// Get the posterior mode 
			    vector<double> ec(p.n, 0);
			    for (int l = 0; l < p.n; l++) {
					ec[l] = (*PsPtr)[l].post_p;
			    }
				int index =0;
				for (int l=0; l <p.n; l++) {
					if(ec[l] > ec[index]) index =l;
				}
				P_bestlevel = (*PsPtr)[index];
	            cout<<"no improvement in 10 levels!\n";
	            break;
	        }
        }else {
        	vector<double> ec(p.n, 0);
			int index = rand_int(ec);
			P_bestlevel = (*PsPtr)[index];
        }	
    }

    //delete newPsPtr; // add by johnmu
    //delete c_originals; // added by johnmu
    delete PsPtr;       //+ by ed520
    delete c_originalsPtr;//+ by ed520
    delete cutregsPtr; //+ by ed520
#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerFinish_cpu( timer.tTime_C_SIS , timer.sTime_C_SIS );
    else
        timer.TimerFinish_cpu(timer.tTime_NC_SIS ,timer.sTime_NC_SIS);
#endif
    return;
}



double Output_distance_2 ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect){

    vector<double>   modewts;
    vector<OnePartition_data> modePs;
    vector<int> modeinds;
    double ess;
    if (pathctcorrect) {
         cout << "P_bestlevel after pathct correction:" << endl;
         for (int par = 0; par < (int) P_bestlevel.size(); par++) {
             P_bestlevel[par].w -= path_log_count(P_bestlevel[par], pathctmap);
         }
    }

    getmax_corr(P_bestlevel, modeinds, modewts, ess);
    modePs.clear();
    for (int lala = 0; lala < (int) modeinds.size(); lala++) {
        modePs.push_back(P_bestlevel[modeinds[lala]]);
    }
 //   vector<double> oneHDs = getvarioushellingerdist(p, smoothneighbordist, density, modePs, modewts);

    vector<double> oneKLs = getvariousKLdist(p, smoothneighbordist, density, modePs, modewts);

 //   cout << " unifHD=" << oneHDs[0] << " HD=" << oneHDs[1] << "\n";
 //   cout << " unifHDsmooth=" << oneHDs[2] << " HDsmooth=" << oneHDs[3] << "\n";

    cout << " unifKL=" << oneKLs[0] << " KL=" << oneKLs[1] << "\n";
    cout << " unifKLsmooth=" << oneKLs[2] << " KLsmooth=" << oneKLs[3] << "\n";
    cout << " unifKLsmooth*2=" << oneKLs[4] << " KLsmooth*2=" << oneKLs[5] << "\n";
    return oneKLs[1];

}

void Output_distance ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect){
        Output_distance_2 (p,P_bestlevel, density, smoothneighbordist, pathctmap,  pathctcorrect);
        return;
}


void Output_distance_mixnormal ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp,  double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect,vector<OnePartition_data> SepPs ){

    vector<double>  modewts;
    vector<int> modeinds;

    vector<OnePartition_data> modePs;
    double ess;

    if (pathctcorrect) {
        cout << "P_bestlevel after pathct correction:" << endl;

        for (int par = 0; par < (int) P_bestlevel.size(); par++) {
            P_bestlevel[par].w -= path_log_count(P_bestlevel[par], pathctmap);
        }
    }
    getmax_corr(P_bestlevel, modeinds, modewts, ess);
    modePs.clear();
    for (int lala = 0; lala < (int) modeinds.size(); lala++) {
        modePs.push_back(P_bestlevel[modeinds[lala]]);
    }

//    vector<double> oneHDs = getvarioushellingerdist_mixnormal(p, smoothneighbordist, density, mixp, modePs, modewts, SepPs);
    vector<double> oneKLs = getvariousKLdist_mixnormal(p, smoothneighbordist, density, mixp, modePs, modewts, SepPs);
 //   cout << " unifHD=" << oneHDs[0] << " HD=" << oneHDs[1] << "\n";
 //   cout << " unifHDsmooth=" << oneHDs[2] << " HDsmooth=" << oneHDs[3] << "\n";

    cout << " unifKL=" << oneKLs[0] << " KL=" << oneKLs[1] << "\n";
    cout << " unifKLsmooth=" << oneKLs[2] << " KLsmooth=" << oneKLs[3] << "\n";

    cout << " smoothneighbordist = " << smoothneighbordist << endl;
//    for (int k = 0; k < (int) optmodeinds.size(); k++) {
//        cout << "wt=" << optmodewts[k] << '\n';
//        print_partition(optmodeinds[k]);
//    }

}

// - by thchiu
/*
// for the plot: BPScore is linear of KL divergence
void SISforlinearity( parameters & p, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp) {

    double* time_array = NULL;
    vector<vector<double> > c_originals;
    c_originals.resize(p.n);
//    vector<double> scores;
    vector<int> cutregs;
    cutregs.resize(p.n);
    vector<double> onelevelscores;
    vector<OnePartition_data> newPs;
    Create(p, newPs);
    vector<OnePartition_data>  EmptySepPs;
    double currlBPscore = - 1;

    for (int i = 0; i < p.levels; i++) {
        cout << "\n\n###i=" << i << endl;
        if ((i > 0) && (p.resampling) && (i % p.steps == 0)) { //(ess<10) &&  with getmax included.
           resample(newPs, p.resampling, c_originals, cutregs);
        }
        //bool exhaust = sample_trees(newPs, c_originals, cutregs, p.maxpercentage, time_array,false);  //commented by thchiu
        sample_trees(newPs, c_originals, cutregs, p.maxpercentage, time_array,false);

        for(int j=0; j<(int)newPs.size(); j++){
        currlBPscore = -lgamma(p.samplesize + (i + 2) / 2.0)- (i + 2) * lgamma(0.5) + lgamma(0.5 * (i + 2)) - 1 * (i + 2);
        currlBPscore = currlBPscore + lBPscore(newPs[j]);

        cout << currlBPscore << "\t";
        cout<< KL_fromf_mixnormal(density, mixp, p.testdata, newPs[j], p.samplesize,0,EmptySepPs)<<endl;
        }

    }
    return;
}
*/

vector<double> Gibbs_FOF_weight(vector<vector<double> >& testdata, vector<OneDensity_data> ensemble, double betaM, double betaP, double splitP, double tree_temperature, double forest_temperature){
	int testsize;
	testsize = testdata.size();
	int dim;
	dim = testdata[0].size();
	//double forest_temperature = 1.0; 
	//double tree_temperature = 50.0; 
	
	vector<double> result(3,0);
	double prior =0.0;
	double tree_prior = 0.0;
	double tree_risk = 0.0;
	
	// keep track of log(density)
	vector<double> est_density(testsize, 0);
	vector<double> valid_count(testsize,0);
	for(int l=0; l<(int)ensemble.size(); l++) {	
		double tree_weights = 0.0; 	
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
			
			// vector<double> shifted_data(testsize, 0);
			vector<double> newtestx(testsize,0);
			for(int i=0; i< testsize; i++) {
				shifted_data[i] =0.5*(testdata[i][d] + ensemble[l].shift[d]);
			}
			
			newtestx = Ftransform(shifted_data, ensemble[l].marginalP[d]);
			
			for(int i=0; i< testsize; i++) {
				// get marginal quantiles for test points based on lth estimate, linxi
				trans_test[i][d] = 0.5*(newtestx[i]+ensemble[l].joint_shift[d]);
			}
			
			// checkings, linxi
			/*
			cerr<<"original test data is: ";
	        for(int i=0; i<10; i++) cerr<<shifted_data[i]<<" "<<testdata[i][d]<<" ";
	        cerr<<endl;

	        cerr<<"transformed test data is: ";
			for(int i=0; i<10; i++) cerr<<newtestx[i]<<" "<<trans_test[i][d]<<" ";
	        cerr<<endl;

	        vector<double> recoveronevartest = inv_Ftransform(newtestx, ensemble[l].marginalP[d]);

			cerr<<"recovered test data is: ";
	        for(int i=0; i<10; i++) cerr<<recoveronevartest[i]<<" "<<2*recoveronevartest[i]-ensemble[l].shift[d]<<" ";
	        cerr<<endl;
			*/
			
			// calculate estimated density based on one tree (product of marginal densities)
			for(int i=0; i < testsize; i++) {
				vector<double> onepiece(1, shifted_data[i]);
				onetreeD[i] *= one_density(ensemble[l].marginalP[d], onepiece, 0.0);
			}
			
			// This is related to the TREE adjustment 
			//double temp_tree_prior = ensemble[l].marginalP[d].prior; // tree prior, after scaling by 2 
			//double temp_tree_risk = ensemble[l].mparginalP[d].BF; // tree risk, after scaling by 2
			
			//int temp_split = ensemble[l].marginalP[d].sregs.size()-1; //The number of splits made to grow the tree
			//temp_tree_prior += (double)temp_size * betaP * log(2.0);// For every split, the depth is counted one level deeper, split prob 0.5 * 2^{-beta *depth}
			//temp_tree_risk -= (double) testsize * log(2.0); // For every individual density, the 
			//tree_weights -= ensemble[l].marginalP[d].post_p;
			tree_weights -= ensemble[l].marginalP[d].BF;
			
		}
		
		for(int i=0; i< testsize; i++) {
			if(dim > 1) {
				// calculate estimated density based on one tree (joint density)
				onetreeD[i]*= one_density(ensemble[l].jointP, trans_test[i],0.0);
				onetreeD[i]*= exp(-(double)(2*dim)*log(2.0));
			
				// update estimated density based on the ensemble
				if(!isnan(onetreeD[i])) {
					est_density[i] += onetreeD[i];
					valid_count[i] +=1.0; 		
				}
			}
			else{
				onetreeD[i]*= exp(-(double)(dim)*log(2.0));
				//est_density[i] += onetreeD[i];
				
				if(!isnan(onetreeD[i])) {
					est_density[i] += onetreeD[i];
					valid_count[i] +=1.0; 		
				}
			}
		}
		
		if(dim >1) {
			//tree_weights -= ensemble[l].jointP.post_p; 
			tree_weights -= ensemble[l].jointP.BF; 
		}
			
		int k=0;
		double tree_risk =0.0;
		for(int i=0; i < testsize; i++) {
			if(onetreeD[i] >0 && !isnan(log(onetreeD[i]))) {
				tree_risk += log(onetreeD[i]);
				k++;
			} else if (onetreeD[i]==0) {
				k++;
			} else{
				cout<<"estimated single tree density is: "<<onetreeD[i]<<endl;
			}
		}
		result[0] = result[0] + tree_weights/tree_temperature;
		
		cout<<"Calculated risk associated with tree "<<l<<" is:"<<tree_risk<<"; Posterior prob is: "<<ensemble[l].marginalP[0].post_p<<"; weight is:"<<ensemble[l].marginalP[0].w<<"; Bayes factor (sequential cal) is:"<<ensemble[l].marginalP[0].BF<<"; Prior is:"<<ensemble[l].marginalP[0].prior<<endl;			
	}
	
	for(int i=0; i < testsize; i++) {
		if(testdata[i][0] <0.00001 || testdata[i][0]>0.99999) {
			est_density[i]=0;
		}else {
			if(valid_count[i] >0 ) {
				est_density[i] *= 1/valid_count[i];
			}else {
				cout<<i<<"\t"<<est_density[i]<<endl;
				est_density[i] = 0; 			
			}
		}		
	}
	
    int k=0;
	double f_weight =0.0;
    for (int i = 0; i < testsize; i++) {
		
        if(est_density[i]>0 && !isnan(log(est_density[i])))  {
            f_weight+= log(est_density[i]);
            k++;
        } else if(est_density[i]==0 ){
        	k++;
        } else{
        	cout<<"estimated density: "<<est_density[i]<<" "<<endl;
        }
		
		if(i < 10) {
			cout<<testdata[i][0]<<" "<<testdata[i][1]<<" "<<"estimated density: "<<est_density[i]<<" "<<endl;
		}
    }
    cout<<"k="<<k<<endl;
    //f_weight /=(double)k;
	result[2] = result[0];
	if(dim ==1) {
		result[0] = result[0] +(1.0/(1-splitP))*f_weight/forest_temperature;	
	}else {
		result[0] = (1.0/(1-splitP))*f_weight/forest_temperature;	
	}  
	result[1] = (1.0/(1-splitP))*f_weight/forest_temperature;
	cout<<"weight associated with forest is:"<<result[0]<<endl;
	cout<<"posterior probability of the forest is:"<<result[1]<<endl;
    return result;

}

void GDF_UQ_sampling(vector<OneForest_data> GDF_posterior, vector<OneForest_data> & UQ_samples, int n_sample, double pseudo_count) {
	int n_forests =(int)GDF_posterior.size();
	vector<double> lweights(n_forests, 0), weights(n_forests, 0);
	double w_lsum =0.0;
	for(int k=0; k<n_forests; k++) {
		lweights[k] = GDF_posterior[k].weight;
		cout<<"The weight associated with forest "<<k<<" is:"<<GDF_posterior[k].weight<<endl;
		w_lsum += lweights[k];
	}
	w_lsum = w_lsum/(double)n_forests;
	cout<<"Average (log) weight is: "<<w_lsum<<endl;
	for(int k=0; k<n_forests; k++) {
		lweights[k] = lweights[k]- w_lsum;
		cout<<"The asjusted weight associated with forest "<<k<<" is:"<<lweights[k]<<endl;
	}
	w_lsum = logsum(lweights);
	cout<<"w_lsum is: "<<w_lsum<<endl;
	double w_sum =0.0;
	for(int k=0; k<n_forests; k++) {
		lweights[k] -= w_lsum;
		weights[k] = exp(lweights[k]);
		w_sum += weights[k];
	}
	cout<<"The summation of normalized weights is: "<<w_sum<<endl;
	
	UQ_samples.clear();
	
	for(int k=0; k< n_sample; k++) {
		// First, sample forest topology
		int sample_index = rand_int(weights);
		UQ_samples.push_back(GDF_posterior[sample_index]);
		
		// Given the forest topology, sample tree-supported densities under the OPT prior
		for(int m=0; m<(int)UQ_samples[k].forestP.size(); m++) {
			// We only consider one-dimensional densities for UQ, therefore, only need to sample one marginalP (first dimension)
			vector<double> sample_density = sample_posterior_density(UQ_samples[k].forestP[m].marginalP[0], pseudo_count);
			if((int)sample_density.size()== (int)UQ_samples[k].forestP[m].marginalP[0].sregs.size()) {
				for(int l=0; l< (int)UQ_samples[k].forestP[m].marginalP[0].sregs.size(); l++) {
					UQ_samples[k].forestP[m].marginalP[0].sregs[l].den = sample_density[l];
				}	
			}else{
				cout<<"An error occurs at the sampling from OPT posterior step: the number of leaf nodes does not match!"<<endl;
				return; 
			}
			
			/* checkings, not needed
			if(k<10 && m<10) {
				cout<<"k= "<<k<<", m= "<<m<<endl;
				vector<double> mle_density = get_all_density(UQ_samples[k].forestP[m].marginalP[0]);
				for(int i=0; i<(int)mle_density.size(); i++) {
					cout<<mle_density[i]<<"\t"<<sample_density[i]<<endl;
					if(i<(int)mle_density.size()-1) {
						cout<<UQ_samples[k].forestP[m].marginalP[0].cut_index[i]<<"\t"<<UQ_samples[k].forestP[m].marginalP[0].left_count[i]<<"\t"<< UQ_samples[k].forestP[m].marginalP[0].right_count[i]<<endl;
					}
				}
				vector<double> prob ={1};
			
				for(int i=0; i<(int)UQ_samples[k].forestP[m].marginalP[0].cut_index.size(); i++) {
					double orig_prob = prob[UQ_samples[k].forestP[m].marginalP[0].cut_index[i]];
					cout<<"i= "<<i<<", cut index: "<<UQ_samples[k].forestP[m].marginalP[0].cut_index[i]<<endl;
					cout<<"orig_prob: "<<orig_prob<<endl;
					int rand_seed = rand_int(0, 65535);
					std::mt19937_64 rng(rand_seed);
					double split_prob = rbeta_gamma_ratio((double)UQ_samples[k].forestP[m].marginalP[0].left_count[i]+ pseudo_count, (double)UQ_samples[k].forestP[m].marginalP[0].right_count[i]+ pseudo_count, rng);
					//double split_prob = rbeta(rand_double(), (double)UQ_samples[k].forestP[m].marginalP[0].left_count[i]+ pseudo_count, (double)UQ_samples[k].forestP[m].marginalP[0].right_count[i]+ pseudo_count);
					cout<<"split_prob: "<<split_prob<<", left: "<<UQ_samples[k].forestP[m].marginalP[0].left_count[i]<<", right: "<<UQ_samples[k].forestP[m].marginalP[0].right_count[i]<<endl;
					prob[UQ_samples[k].forestP[m].marginalP[0].cut_index[i]] = orig_prob * split_prob;
					prob.push_back (orig_prob * (1-split_prob));
				}
			}
			*/			
		}		
	}
}

void GDF_UQ_path_evaluation(vector<OneForest_data> & UQ_samples, vector<vector<double>> test_points, string out_dir) {
	int n_sample = (int)UQ_samples.size();
	int dim = (int) test_points[0].size();
	int testsize = (int)test_points.size();
	
	string cline = "mkdir -p "+out_dir;
	const int dir= system(cline.c_str());
	if (dir< 0) return;
	
	for(int k=0; k< n_sample; k++) {
		// output for the kth forest
		string ofilenamerep = out_dir + "/UQforest_"+std::to_string(k)+"/";
		string clinerep = "mkdir -p "+ofilenamerep;
		const int dir= system(clinerep.c_str());
		if (dir< 0) return;
		
		string outdenfile = ofilenamerep+"forest_density.txt";
		ofstream outputrep(outdenfile.c_str());
		
		vector<double> oneforestD(testsize,0);
		for(int l=0; l<(int) UQ_samples[k].forestP.size(); l++) {		
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
				cerr<<"shift along dimension "<<d<<" is: "<<UQ_samples[k].forestP[l].shift[d]<<endl;
			
				vector<double> newtestx(testsize,0);
				for(int i=0; i< testsize; i++) {
					shifted_data[i] =0.5*(test_points[i][d] + UQ_samples[k].forestP[l].shift[d]);
				}
			
				/* Not needed in the 1-D case
				newtestx = Ftransform(shifted_data, UQ_samples[k].forestP[l].marginalP[d]);
			
				for(int i=0; i< testsize; i++) {
					// get marginal quantiles for test points based on lth estimate, linxi
					trans_test[i][d] = 0.5*(newtestx[i]+UQ_samples[k].forestP[l].joint_shift[d]);
				}
				*/
			
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
					onetreeD[i] *= one_density_UQ(UQ_samples[k].forestP[l].marginalP[d], onepiece);
					/*
					if(i<10 && l<10) {
						cout<<"The shifted data is: "<<shifted_data[i]<<"; The estimated density is: "<<one_density_UQ(UQ_samples[k].forestP[l].marginalP[d], onepiece)<<endl;
					}
					*/
				}
			
			}
		
			for(int i=0; i< testsize; i++) {
				if(dim > 1) {
					/* Not needed in the 1-D case
					// calculate estimated density based on one tree (joint density)
					onetreeD[i]*= one_density(UQ_samples[k].forestP[l].jointP, trans_test[i], 0.0);
					onetreeD[i]*= exp(-(double)(2*dim)*log(2.0));
			
					// update estimated density based on the ensemble
					oneforestD[i] += onetreeD[i];
					*/
				}
				else{
					onetreeD[i]*= exp(-(double)(dim)*log(2.0));
					oneforestD[i] += onetreeD[i];
				}
			}			
		}
		for(int i=0; i<testsize; i++) {
			oneforestD[i] *= 1/(double)UQ_samples[k].forestP.size();
		}	
		// output the evaluated denstifies at test points
		for(int i=0; i<testsize; i++) {
			for(int d=0; d <dim; d++) {
				outputrep<< test_points[i][d]<<"\t";
			}
			outputrep << oneforestD[i]<<endl;
		}
		outputrep.close();
	}
}
