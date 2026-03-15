/*
start from \Omega, N trees
1. update trees to add one more cut (based on inexact OPT)
2. calculate w_t = w_t-1*gamma(0.5+n_k^1)gamma(0.5+n_k^2)/gamma(0.5+n_k)(\mu(A_k^1)\mu(A_k^2)/\mu(A_k))^\beta|A_k|^n_k/||||/\alpha_k\lambda_j
3. if necessary, do resampling
*/




#include "sampling.h"
#include "output.h"
#include "count.h"  //+ by thchiu
#include "timer_ed520.h" //+ by ed520 
extern TIMER timer; //+ by ed520

// consider the discrete case: for any region, it will never cut the dimension in which the value of the sample points all agree. If there is no region no dimension to be cut, 'c' will be empty.
// P has t regions: P.sregs.size()=t, then c_original before excuting this function is of length dim*(t-1), cutreg=0,1,...,t-2
// c_original are the original Bayes factors, comparable across different levels, linxi
// after excuting this function, c_original[t*dim +d] is the Bayes factor for further splitting region i along dimension d
// c are Bayes factors for regions, c[regionID*dim+d] is the BF for splitting region indexed by regionID along dimension d
// c is orginalized as (region1, dim1) ... (region1, dim d), (region2, dim 1), ... region(nRegion, dim d)
// c are weights after scaling, according to "smooth" parameter, used for splitting to avoid numerical issue, linxi
//void get_info_from_onetree(OnePartition_data &P, vector<double> & c_original,vector<double> &c,int cutreg, double maxpercentage, double &smooth, double *time_array, bool discrete) { //- by thchiu
void get_info_from_onetree(OnePartition_data &P, vector<double> & c_original,vector<double> &c,int cutreg, double maxpercentage, double &smooth, double *time_array, bool discrete, parameters & para, int Cdepth, string method){ //+ by thchiu
#ifdef PROFILE
    if(para.dim == 1)
        timer.TimerStart_cpu(timer.sTime_get_info_C);
    else
        timer.TimerStart_cpu(timer.sTime_get_info_NC);
#endif
	/*
    //double beta=0.5;                //pseudo-ct for dirichlet. //- by thchiu
    int dim = P.sregs[0].reg_code.size();

    int Tmp_count[dim];     //+ by thchiu
	int Tmp_count_vali[dim];
    bool Tmp_unique[dim];   //+ by thchiu

    if ((int) c_original.size() == 0) {
        for (int d=0; d<dim; d++)   c_original.push_back(0);

        Count(P.sregs[0], Tmp_count, Tmp_count_vali, Tmp_unique, dim, discrete, para);    //+ by thchiu

        Calculate_Weight(c_original, Tmp_unique, dim, Tmp_count, P.sregs[0].num, Tmp_count_vali, P.sregs[0].valiNum, 0);  //+ by thchiu
    }


    else {

        if ((int) c_original.size() != ((int) P.sregs.size() - 1) * dim) {
            cout<<" c_original.size() ="<< c_original.size() <<endl;
            cout<<c_original[0]<<endl;
            cout<<"P.sregs.size()="<<P.sregs.size()<<endl;
            cout << "get_info wrong c_original!\n";
            exit(5);
        }

        for (int d=0; d<dim; d++)   c_original.push_back(0);
        //shrink_reg_data_short & sreg = P.sregs[cutreg];   //- by thchiu

        //Below 4 lines are added by thchiu
        Count(P.sregs[cutreg], Tmp_count, Tmp_count_vali, Tmp_unique, dim, discrete, para);
        Calculate_Weight(c_original, Tmp_unique, dim, Tmp_count, P.sregs[cutreg].num, Tmp_count_vali, P.sregs[cutreg].valiNum, cutreg);
        Count(P.sregs[(int) P.sregs.size() - 1], Tmp_count, Tmp_count_vali, Tmp_unique, dim, discrete, para);
        Calculate_Weight(c_original, Tmp_unique, dim, Tmp_count, P.sregs[(int) P.sregs.size() - 1].num, Tmp_count_vali, P.sregs[(int) P.sregs.size() - 1].valiNum, (int) P.sregs.size() - 1);
    }
	
	int nRegion =0; // number of coupled regions in current tree P
	
	for(int i=0; i<(int)P.sregs.size(); i++) {
		int tempR = P.sregs[i].regionID;
		if(tempR > nRegion) 
			nRegion = tempR;
	}
	nRegion +=1; //as the index of region start from 0, the number of regions should be nRegion after for loop +1
	
	// c is orginalized as (region1, dim1) ... (region1, dim d), (region2, dim 1), ... region(nRegion, dim d)
	
	c.resize((int)nRegion*dim);
	for(int i=0; i< (int)c.size(); i++) {
		c[i]=0;
	}
	for(int i=0; i<(int)P.sregs.size(); i++) {
		int tempID = P.sregs[i].regionID;
		for( int d=0; d< dim; d++) {
			c[tempID*dim+d] += c_original [i*dim +d]; // will add log BF for all nodes within the same region together
		}
	}
	
	// This part of code is for checking elements of c
    double cmax= c[0];
    int maxind= 0;

     for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] > cmax ){
            cmax = c[i];
            maxind=i;
        }
    }

    if (cmax <1 -MAXNUM) {
        c.clear();
        cout << "no regions no dimensions to be split!!" << endl;
      return;
    }
    double cmax2 = -2 * MAXNUM; // cmax2 is the second maximum of c
    for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] > cmax2 && c[i] < cmax) cmax2 = c[i];
    }
    if (cmax2 <1 -MAXNUM) return;

   for(int i = 0; i < (int) c.size(); i++) {
       if(isnan(c[i])) {
                     cout<<"i="<<i<<"c[i]=NAN at a"<<endl;
       }
    }
	
	// smoothing to make sure weights are not too small, linxi
    if (maxpercentage < 0) smooth = 20;
    else smooth = (cmax - cmax2) / log(maxpercentage / (1 - maxpercentage)); // smooth=1: no smooth
    if(isinf(smooth)) cout<<"smooth==inf! cmax="<<cmax<<" cmax2="<<cmax2<<endl;
    if (smooth > 1) {
        for (int i = 0; i < (int) c.size(); i++) {
            if( c[i] > 1-MAXNUM )           c[i] = c[i] / smooth;
        }
    }
     for(int i = 0; i < (int) c.size(); i++) {
       if(isnan(c[i])) {
                     cout<<"i="<<i<<"c[i]=NAN at b"<<endl;
                     cout<<"smooth="<<smooth<<endl;
       }
    }
	*/
    //double beta=0.5;                //pseudo-ct for dirichlet. //- by thchiu
    int dim = P.sregs[0].reg_code.size();

    int Tmp_count[dim];     //+ by thchiu
	int Tmp_count_vali[dim];
    bool Tmp_unique[dim];   //+ by thchiu

    if ((int) c_original.size() == 0) {
        for (int d=0; d<dim; d++)   c_original.push_back(0);
        //vector<usint_mask> new_reg_code;  //- by thchiu
        //int i = 0;    //- by thchiu
        //shrink_reg_data_short sreg = P.sregs[i]; // if the ith region is updated  // - by thchiu

/*{{{*/
        /*for (int d = 0; d < dim; d++) { // if the ith region is splitted according to the dth dimension
            new_reg_code = sreg.reg_code;
            beta=0.5;
            if ((double) sreg.num / 200.0 < 0.5) beta = max(0.1, (double) sreg.num / 200.0); // heuristic pseudo ct.
            c_original[dim * i + d] = 0;
            for (int j = 0; j < 2; j++) {   //problem10: 1 loop is enough, thchiu
                new_reg_code[d].x = (sreg.reg_code[d].x << 1) + j; // generate a child's regcode
                new_reg_code[d].mask = (sreg.reg_code[d].mask << 1) + 1;
                bool unique = false;
                shrink_reg_data_short sr;
                PtsInReg_short(sr, new_reg_code, sreg.regdata, sreg.num, d, unique, false, discrete);

                if (unique) {// || new_reg_code[d].mask >= 16) {  second part : for letter data(cont)
                    c_original[dim * i + d] = -MAXNUM;
                    break;
                }
                c_original[dim * i + d] += lgamma(beta + sr.num);
            }
            if (c_original[dim * i + d] > 1 - MAXNUM) c_original[dim * i + d] -= lgamma(beta + sreg.num);
            if (c_original[dim * i + d] > 1 - MAXNUM) c_original[dim * i + d] += sreg.num * log(2.0);
        }*//*}}}*/

        Count(P.sregs[0], Tmp_count, Tmp_count_vali, Tmp_unique, dim, discrete, para);    //+ by thchiu
		
		if(method =="gibbs") {
			Calculate_Weight_Gibbs(c_original, Tmp_unique, dim, Tmp_count, P.sregs[0].num, Tmp_count_vali, P.sregs[0].valiNum, 0);  //+ by thchiu		
		}else if(method == "opt") {
			Calculate_Weight_OPT(c_original, Tmp_unique, dim, Tmp_count, P.sregs[0].num, Tmp_count_vali, P.sregs[0].valiNum, 0); 
		}else if(method == "RF") {
			Calculate_Weight_RF(c_original, Tmp_unique, dim, Tmp_count, P.sregs[0].num, Tmp_count_vali, P.sregs[0].valiNum, 0); 
		}else {
			cout<<"The method for tree growth is not specified!"<<endl;
		}
        
    }


    else {

        if ((int) c_original.size() != ((int) P.sregs.size() - 1) * dim) {
            cout<<" c_original.size() ="<< c_original.size() <<endl;
            cout<<c_original[0]<<endl;
            cout<<"P.sregs.size()="<<P.sregs.size()<<endl;
            cout << "get_info wrong c_original!\n";
            exit(5);
        }

        for (int d=0; d<dim; d++)   c_original.push_back(0);
        //shrink_reg_data_short & sreg = P.sregs[cutreg];   //- by thchiu

        //Below 4 lines are added by thchiu
        Count(P.sregs[cutreg], Tmp_count, Tmp_count_vali, Tmp_unique, dim, discrete, para);
		if(method == "gibbs") {
			Calculate_Weight_Gibbs(c_original, Tmp_unique, dim, Tmp_count, P.sregs[cutreg].num, Tmp_count_vali, P.sregs[cutreg].valiNum, cutreg);	
		}else if(method == "opt") {
			Calculate_Weight_OPT(c_original, Tmp_unique, dim, Tmp_count, P.sregs[cutreg].num, Tmp_count_vali, P.sregs[cutreg].valiNum, cutreg);
		}else if(method == "RF") {
			Calculate_Weight_RF(c_original, Tmp_unique, dim, Tmp_count, P.sregs[cutreg].num, Tmp_count_vali, P.sregs[cutreg].valiNum, cutreg);
			// for RF, stop partitioning when the depth of the region reaches the limit
			int depth = get_depth(P.sregs[cutreg].reg_code);
			if(depth == Cdepth) {
				int RegionID = cutreg;
			    for(int i=0; i < dim; ++i){
					// uniform weight, randomly choose a direction to split. 
			        c_original[dim*RegionID+i] = -MAXNUM;
			    }
			}
		}else {
			cout<<"The method for tree growth is not specified!"<<endl;
		}
		
        Count(P.sregs[(int) P.sregs.size() - 1], Tmp_count, Tmp_count_vali, Tmp_unique, dim, discrete, para);
		if(method == "gibbs") {
			Calculate_Weight_Gibbs(c_original, Tmp_unique, dim, Tmp_count, P.sregs[(int) P.sregs.size() - 1].num, Tmp_count_vali, P.sregs[(int) P.sregs.size() - 1].valiNum, (int) P.sregs.size() - 1);
		}else if(method == "opt") {
			Calculate_Weight_OPT(c_original, Tmp_unique, dim, Tmp_count, P.sregs[(int) P.sregs.size() - 1].num, Tmp_count_vali, P.sregs[(int) P.sregs.size() - 1].valiNum, (int) P.sregs.size() - 1);
		}else if(method == "RF") {
			Calculate_Weight_RF(c_original, Tmp_unique, dim, Tmp_count, P.sregs[(int) P.sregs.size() - 1].num, Tmp_count_vali, P.sregs[(int) P.sregs.size() - 1].valiNum, (int) P.sregs.size() - 1);
			// for RF, stop partitioning when the depth of the region reaches the limit
			int depth = get_depth(P.sregs[(int) P.sregs.size() - 1].reg_code);
			if(depth == Cdepth) {
				int RegionID = (int) P.sregs.size() - 1;
			    for(int i=0; i < dim; ++i){
					// for regions reaching the depth limit, the probability of making a further split is close to 0
			        c_original[dim*RegionID+i] = -MAXNUM;
			    }
			}
		}else {
			cout<<"The method for tree growth is not specified!"<<endl;
		}
        

        /*{{{*/
        /*vector<usint_mask> new_reg_code;      //- by thchiu
        vector<int> updateindex;
        updateindex.push_back(cutreg);
        updateindex.push_back((int) P.sregs.size() - 1);

        for (int kk = 0; kk < 2; kk++) {
            int i = updateindex[kk];


            //shrink_reg_data_short sreg(P.sregs[i]); // if the ith region is updated   //- by thchiu
            shrink_reg_data_short & sreg = P.sregs[i]; // + by thchiu

            for (int d = 0; d < dim; d++) { // if the ith region is splitted according to the dth dimension
                new_reg_code = sreg.reg_code;
                beta=0.5;
                if ((double) sreg.num / 200.0 < 0.5) beta = max(0.1, (double) sreg.num / 200.0); // heuristic pseudo ct.
                c_original[dim * i + d] = 0;
                for (int j = 0; j < 2; j++) {
                    new_reg_code[d].x = (sreg.reg_code[d].x << 1) + j; // generate a child's regcode
                    new_reg_code[d].mask = (sreg.reg_code[d].mask << 1) + 1;
                    bool unique = false;
                    shrink_reg_data_short  sr;

                    PtsInReg_short(sr, new_reg_code, sreg.regdata, sreg.num, d, unique, false, discrete);

                    if (unique) {// || new_reg_code[d].mask >= 16) {  second part : for letter data(cont)
                        c_original[dim * i + d] = -MAXNUM;
                        break;
                    }
                    c_original[dim * i + d] += lgamma(beta + sr.num);
                }
                if (c_original[dim * i + d] > 1 - MAXNUM) c_original[dim * i + d] -= lgamma(beta + sreg.num);
                if (c_original[dim * i + d] > 1 - MAXNUM) c_original[dim * i + d] += sreg.num * log(2.0);
            }

        }*//*}}}*/

    }
    c = c_original;
    double cmax= c[0];
    int maxind= 0;

     for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] > cmax ){
            cmax = c[i];
            maxind=i;
        }
    }

    if (cmax <1 -MAXNUM) {
        c.clear();
        cout << "no regions no dimensions to be split!!" << endl;
      return;
    }
	for (int i=0; i< (int)c.size(); i++) {
		c[i] -= cmax;
	}
	cmax =0;
	
	// smoothing to make sure weights are not too small, linxi
	// only for gibbs or opt
	if(method != "RF") {
	    double cmax2 = - MAXNUM; // cmax2 is the second maximum of c
	    for (int i = 0; i < (int) c.size(); i++) {
	        if (c[i] > cmax2 && c[i] < cmax) cmax2 = c[i];
	    }
	    if (cmax2 <1 -MAXNUM) return;

	   for(int i = 0; i < (int) c.size(); i++) {
	       if(isnan(c[i])) {
	                     cout<<"i="<<i<<"c[i]=NAN at a"<<endl;
	       }
	    }
	
	    if (maxpercentage < 0) smooth = 20;
	    else smooth = (cmax - cmax2) / log(maxpercentage / (1 - maxpercentage)); // smooth=1: no smooth
	    if(isinf(smooth)) cout<<"smooth==inf! cmax="<<cmax<<" cmax2="<<cmax2<<endl;
	    if (smooth > 1) {
	        for (int i = 0; i < (int) c.size(); i++) {
	            if( c[i] > 1-MAXNUM )           c[i] = c[i] / smooth;
	        }
	    }
		for(int i = 0; i < (int) c.size(); i++) {
			if(isnan(c[i])) {
				cout<<"i="<<i<<"c[i]=NAN at b"<<endl;
				cout<<"smooth="<<smooth<<endl;
	       }
	   }	
	}
#ifdef PROFILE
     if(para.dim == 1)
         timer.TimerFinish_cpu(timer.tTime_get_info_C , timer.sTime_get_info_C);
     else
         timer.TimerFinish_cpu(timer.tTime_get_info_NC , timer.sTime_get_info_NC);
#endif
}

//void rand_choose_one_nextleveltree(OnePartition_data &newP, OnePartition_data &P, vector<double> &c, int &cutreg, double smooth, bool discrete) { //- by thchiu
void rand_choose_one_nextleveltree(OnePartition_data &P, vector<double> & c_original, vector<double> &c, int &cutreg, double smooth, bool discrete, parameters& para, double beta, int Cdepth, string method, double tree_temperature, bool UQ_flag) {// + by thchiu
#ifdef PROFILE
    if(para.dim == 1)
        timer.TimerStart_cpu(timer.sTime_rand_choose_C);
    else
        timer.TimerStart_cpu(timer.sTime_rand_choose_NC);
#endif
	/*
    int dim = P.sregs[0].reg_code.size();
    vector<double> ec(dim * P.sregs.size(), 0);
    double cmax = max(c);
    for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] < 1 - MAXNUM) {
            ec[i] = 0;
        } else ec[i] = exp(c[i] - cmax + 5);
        if (isnan(c[i])) cout << "is.nan(ec[" << i << "]), c[i]=" << c[i] << endl;
    }

    int cut = rand_int(ec);
    int cutvar = cut % dim;
    cutreg = (cut - cutvar) / dim;
	
	int cutdepth;
	vector<int> cutregID; // multiple nodes can be cut at the same time, this is the indices for regions to be cut
	
	cutregID.clear();
	for(int i=0; i<(int)P.sregs.size(); i++){
		if(P.sregs[i].regionID == cutreg) {
			cutdepth = get_depth(P.sregs[i].reg_code); // cutdepth is the depth of the selected region before making the split
			cutregID.push_back(i);
		}
	}
	
	int nRegion = (int)c.size()/dim; // total number of current regions
	
	for(int i=0; i<(int)cutregID.size(); i++) {
	    //newP = P;       //problem8: use P is OK.  - by thchiu
	    //shrink_reg_data_short sreg = P.sregs[cutreg];   //problem8, - by thchiu
	    shrink_reg_data_short & sreg = P.sregs[cutregID[i]]; //reference, + by thchiu
	    //vector<usint_mask> new_reg_code = sreg.reg_code;//problem8, - by thchiu

	    //Begin + by thchiu
	    double Tmp_upbound;
	    sreg.reg_code[cutvar].x = (sreg.reg_code[cutvar].x << 1);
	    sreg.reg_code[cutvar].mask = (sreg.reg_code[cutvar].mask << 1) + 1;
	    pair<double, double> ranges = convert_ranges(sreg.reg_code[cutvar]);
	    Tmp_upbound=ranges.second;

	    shrink_reg_data_short sr;

	    int data_head = sreg.region_start;
	    int data_tail = data_head+sreg.num;
	    int Tmp_Count_train=0;
	    for(;data_head != data_tail;){
	        if(para.data1D[para.pt_start[data_head]+cutvar] > Tmp_upbound){
	            --data_tail;
	            int tmp=para.pt_start[data_head];
	            para.pt_start[data_head]=para.pt_start[data_tail];
	            para.pt_start[data_tail]=tmp;
	            //Tmp_Count++;
	        }
	        else{
	            Tmp_Count_train++;
	            ++data_head;
	        }
	    }

	    if(data_head < sreg.region_start + sreg.num)
	        if(para.data1D[para.pt_start[data_head]+cutvar] <= Tmp_upbound) // + by ed520
	            Tmp_Count_train++;
	    //sr.region_start = sreg.region_start + sreg.regdata.size(); // + by ed520
	    sr.region_start= sreg.region_start + Tmp_Count_train;
	
	    data_head = sreg.region_start_vali;
	    data_tail = data_head+sreg.valiNum;
	    int Tmp_Count_vali=0;
	    for(;data_head != data_tail;){
	        if(para.data1D_vali[para.pt_start_vali[data_head]+cutvar] > Tmp_upbound){
	            --data_tail;
	            int tmp=para.pt_start_vali[data_head];
	            para.pt_start_vali[data_head]=para.pt_start_vali[data_tail];
	            para.pt_start_vali[data_tail]=tmp;
	            //Tmp_Count++;
	        }
	        else{
	            Tmp_Count_vali++;
	            ++data_head;
	        }
	    }

	    if(data_head < sreg.region_start_vali + sreg.valiNum)
	        if(para.data1D_vali[para.pt_start_vali[data_head]+cutvar] <= Tmp_upbound) // + by ed520
	            Tmp_Count_vali++;
	    //sr.region_start = sreg.region_start + sreg.regdata.size(); // + by ed520
	    sr.region_start_vali= sreg.region_start_vali+ Tmp_Count_vali;

	    //sr.num = sr.regdata.size();
	    sr.num = sreg.num - Tmp_Count_train;
		sr.valiNum = sreg.valiNum -Tmp_Count_vali;
	    sr.reg_code= sreg.reg_code;
	    sr.reg_code[cutvar].x ++;
	
	

	    //sreg.num = sreg.regdata.size();
	    sreg.num = Tmp_Count_train;
		sreg.valiNum = Tmp_Count_vali;
		
		if(cutdepth >= Cdepth) {
			sr.regionID = sreg.regionID;  // This is the case that we cut a high-level region, coupling is needed. 
		}else{
			sr.regionID = nRegion; // This is the case NO coupling is needed, the new node is also a new region. Note, index of the region start from 0. 
		}
	    P.sregs.push_back(sr);
	    //End + by thchiu		
	}
 
	double Bfactor =0;
	for(int i=0; i<(int)cutregID.size(); i++) {
		Bfactor += c_original[cutregID[i]*dim + cutvar];
	}
    //double Bfactor = c_original[cut];
	
	//int depth = get_depth(sr.reg_code);
	//cout<<"The depth of region "<<cutreg<<" is: "<<depth-1<<endl;
	
	
	//double priorR = - log((double) dim);
	//// opt stopping rule: tau = 0.5* 2^{- 1.0 * depth}, linxi
	//double stop_rule =0.0;
	//stop_rule = log(0.5) - beta* ((double) (depth-1) )*log(2.0);
	//priorR -= log(1- exp(stop_rule));
	//priorR += stop_rule;
	//stop_rule -= beta*log(2.0);
	//priorR += 2* log(1- exp(stop_rule));
	
	double priorR = - log((double)dim);
	// opt stopping rule: tau = 0.5* 2^{- beta * depth}, linxi
	double stop_rule =log(0.5) - beta* ((double) cutdepth )*log(2.0);
	priorR -= log(1- exp(stop_rule));
	priorR += stop_rule;
	stop_rule -= beta*log(2.0);
	if(cutdepth >= Cdepth) {
		priorR += log(1- exp(stop_rule)); // as regions are coupled, there is no "2*".
	}else {
		priorR += 2*log(1- exp(stop_rule));
	}
	
	// update the weights m
    //if (smooth > 1) newP.w = P.w + (smooth - 1) * c[cut] + logsum(c); // + logsum(c) ;  //problem9: after fixing problem8, this can be P.w = P.w +...., - by thchiu
    // else newP.w = P.w + logsum(c);  //problem9, thchiu
	// "if" part is the weight after smoothing, see derivation, linxi
	P.w = P.w + priorR + Bfactor - c[cut]+ logsum(c);
	P.post_p = P.post_p  + priorR+ Bfactor;
	P.BF = P.BF + Bfactor;
    //if (smooth > 1) P.w = P.w + (smooth - 1) * c[cut] + logsum(c); // + logsum(c); + by thchiu
    //else P.w = P.w + logsum(c); //+ by thchiu
    //return (newP);
	*/
	
    int dim = P.sregs[0].reg_code.size();
    vector<double> ec(dim * P.sregs.size(), 0);
    double cmax = max(c);
    for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] < 1 - MAXNUM) {
            ec[i] = 0;
        } else ec[i] = exp(c[i] - cmax + 5);
        if (isnan(c[i])) cout << "is.nan(ec[" << i << "]), c[i]=" << c[i] << endl;
    }

    int cut = rand_int(ec);
    int cutvar = cut % dim;
    cutreg = (cut - cutvar) / dim;

    //newP = P;       //problem8: use P is OK.  - by thchiu
    //shrink_reg_data_short sreg = P.sregs[cutreg];   //problem8, - by thchiu
    shrink_reg_data_short & sreg = P.sregs[cutreg]; //reference, + by thchiu
    //vector<usint_mask> new_reg_code = sreg.reg_code;//problem8, - by thchiu

/*{{{*/
    //Below is commented by thchiu
    /*for (int j = 0; j < 2; j++) {   //problem7: count 2 time, 1 time is enough, thchiu
        new_reg_code[cutvar].x = (sreg.reg_code[cutvar].x << 1) + j; // generate a child's regcode
        new_reg_code[cutvar].mask = (sreg.reg_code[cutvar].mask << 1) + 1;
        //     if( new_reg_code[cutvar].mask >=16){
        //         cout<<" new_reg_code[cutvar].mask= "<< new_reg_code[cutvar].mask<<"c="<<c[cutreg*dim+cutvar]<<endl;
        //   }
        bool unique = false;
        shrink_reg_data_short sr;
        PtsInReg_short(sr, new_reg_code, sreg.regdata, sreg.num, cutvar, unique, true, discrete); ////

        //if (j == 0) newP.sregs[cutreg] = sr;  //- by thchiu
        //else newP.sregs.push_back(sr);
        if (j == 0) P.sregs[cutreg] = sr;  //+ by thchiu
        else P.sregs.push_back(sr); //+ by thchiu

    }*//*}}}*/

    //Begin + by thchiu
    double Tmp_upbound;
    sreg.reg_code[cutvar].x = (sreg.reg_code[cutvar].x << 1);
    sreg.reg_code[cutvar].mask = (sreg.reg_code[cutvar].mask << 1) + 1;
    pair<double, double> ranges = convert_ranges(sreg.reg_code[cutvar]);
    Tmp_upbound=ranges.second;

    shrink_reg_data_short sr;
    /*for(unsigned i=0; i<sreg.regdata.size(); i++){
        if(sreg.regdata[i][cutvar] > Tmp_upbound){
            sr.regdata.push_back(sreg.regdata[i]);
            sreg.regdata[i] = sreg.regdata[sreg.regdata.size()-1];
            sreg.regdata.pop_back();
            i--;
        }
    }*/

    int data_head = sreg.region_start;
    int data_tail = data_head+sreg.num;
    int Tmp_Count_train=0;
    for(;data_head != data_tail;){
        if(para.data1D[para.pt_start[data_head]+cutvar] > Tmp_upbound){
            --data_tail;
            int tmp=para.pt_start[data_head];
            para.pt_start[data_head]=para.pt_start[data_tail];
            para.pt_start[data_tail]=tmp;
            //Tmp_Count++;
        }
        else{
            Tmp_Count_train++;
            ++data_head;
        }
    }

    if(data_head < sreg.region_start + sreg.num)
        if(para.data1D[para.pt_start[data_head]+cutvar] <= Tmp_upbound) // + by ed520
            Tmp_Count_train++;
    //sr.region_start = sreg.region_start + sreg.regdata.size(); // + by ed520
    sr.region_start= sreg.region_start + Tmp_Count_train;
	
    data_head = sreg.region_start_vali;
    data_tail = data_head+sreg.valiNum;
    int Tmp_Count_vali=0;
    for(;data_head != data_tail;){
        if(para.data1D_vali[para.pt_start_vali[data_head]+cutvar] > Tmp_upbound){
            --data_tail;
            int tmp=para.pt_start_vali[data_head];
            para.pt_start_vali[data_head]=para.pt_start_vali[data_tail];
            para.pt_start_vali[data_tail]=tmp;
            //Tmp_Count++;
        }
        else{
            Tmp_Count_vali++;
            ++data_head;
        }
    }

    if(data_head < sreg.region_start_vali + sreg.valiNum)
        if(para.data1D_vali[para.pt_start_vali[data_head]+cutvar] <= Tmp_upbound) // + by ed520
            Tmp_Count_vali++;
    //sr.region_start = sreg.region_start + sreg.regdata.size(); // + by ed520
    sr.region_start_vali= sreg.region_start_vali+ Tmp_Count_vali;
	
	
	int Tmp_Count_infer=0;
	if(para.dim ==1 && UQ_flag) {
	    data_head = sreg.region_start_infer;
	    data_tail = data_head+sreg.inferNum;
	    for(;data_head != data_tail;){
	        if(para.data1D_infer[para.pt_start_infer[data_head]+cutvar] > Tmp_upbound){
	            --data_tail;
	            int tmp=para.pt_start_infer[data_head];
	            para.pt_start_infer[data_head]=para.pt_start_infer[data_tail];
	            para.pt_start_infer[data_tail]=tmp;
	            //Tmp_Count++;
	        }
	        else{
	            Tmp_Count_infer++;
	            ++data_head;
	        }
	    }

	    if(data_head < sreg.region_start_infer + sreg.inferNum)
	        if(para.data1D_infer[para.pt_start_infer[data_head]+cutvar] <= Tmp_upbound) // + by ed520
	            Tmp_Count_infer++;
	    //sr.region_start = sreg.region_start + sreg.regdata.size(); // + by ed520
	    sr.region_start_infer= sreg.region_start_infer+ Tmp_Count_infer;
	}

    //sr.num = sr.regdata.size();
    sr.num = sreg.num - Tmp_Count_train;
	sr.valiNum = sreg.valiNum -Tmp_Count_vali;
	if(para.dim == 1 && UQ_flag) {
		sr.inferNum = sreg.inferNum -Tmp_Count_infer;
	}
    sr.reg_code= sreg.reg_code;
    sr.reg_code[cutvar].x ++;
	
	

    //sreg.num = sreg.regdata.size();
    sreg.num = Tmp_Count_train;
	sreg.valiNum = Tmp_Count_vali;
	if(para.dim ==1 && UQ_flag) {
		sreg.inferNum = Tmp_Count_infer;
	}
    P.sregs.push_back(sr);
    //End + by thchiu
	
	P.cut_index.push_back(cutreg);
	if(UQ_flag) {
		P.left_count.push_back(Tmp_Count_infer);
		P.right_count.push_back(sr.inferNum);
	}	
	
	if(method == "opt") {
	    double Bfactor = c_original[cut];
	
		int depth = get_depth(sr.reg_code);
		//cout<<"The depth of region "<<cutreg<<" is: "<<depth-1<<endl;
	
		double priorR = - log((double) dim);
		// opt stopping rule: tau = 0.5* 2^{- 1.0 * depth}, linxi
		double stop_rule =0.0;
		stop_rule = log(0.5) - beta * ((double) (depth-1) )*log(2.0);
		priorR -= log(1- exp(stop_rule));
		priorR += stop_rule;
		stop_rule -= beta *log(2.0);
		priorR += 2* log(1- exp(stop_rule));
	
		// update the weights m
	    /*if (smooth > 1) newP.w = P.w + (smooth - 1) * c[cut] + logsum(c); // + logsum(c) ;  //problem9: after fixing problem8, this can be P.w = P.w +...., - by thchiu
	    else newP.w = P.w + logsum(c);  //problem9, thchiu*/
		// "if" part is the weight after smoothing, see derivation, linxi
		P.w = P.w + priorR + Bfactor - c[cut]+ logsum(c);
		P.post_p = P.post_p  + priorR+ Bfactor;
		P.BF = P.BF + Bfactor;
		P.prior += priorR;
	    //if (smooth > 1) P.w = P.w + (smooth - 1) * c[cut] + logsum(c); // + logsum(c); + by thchiu
	    //else P.w = P.w + logsum(c); //+ by thchiu
	    //return (newP);
		if( depth > P.depth) P.depth +=1; 
	}else if(method == "gibbs") {
		//double tree_temperature = 50.0; 
	    double Bfactor = c_original[cut];
	
		int depth = get_depth(sr.reg_code);
		//cout<<"The depth of region "<<cutreg<<" is: "<<depth-1<<endl;
	
		double priorR = - log((double) dim);
		// opt stopping rule: tau =  2^{- beta (1.0 as default)* (depth+1)}, linxi
		/* the stop rule that decays exponentially fast as the depth increases
		double stop_rule =0.0;
		stop_rule =  - beta * ((double) (depth) )*log(2.0);
		priorR -= log(1- exp(stop_rule));
		priorR += stop_rule;
		stop_rule -= beta *log(2.0);
		priorR += 2* log(1- exp(stop_rule));
		*/
		
		// constant stop/split probability 
		// split probability: beta
		priorR -= log(1- beta);
		priorR += log (beta);
		//stop_rule -= beta *log(2.0);
		priorR += 2* log(1- beta);
		
		// This is BSP prior 
		priorR = 0.0;
		priorR -= 1.0/ tree_temperature; 
	
		// update the weights m
	    /*if (smooth > 1) newP.w = P.w + (smooth - 1) * c[cut] + logsum(c); // + logsum(c) ;  //problem9: after fixing problem8, this can be P.w = P.w +...., - by thchiu
	    else newP.w = P.w + logsum(c);  //problem9, thchiu*/
		// "if" part is the weight after smoothing, see derivation, linxi
		P.w = P.w + priorR + Bfactor/tree_temperature- c[cut]+ logsum(c);
		P.post_p = P.post_p  + priorR+ Bfactor/tree_temperature;
		P.BF = P.BF + Bfactor;
		P.prior += priorR;
	    //if (smooth > 1) P.w = P.w + (smooth - 1) * c[cut] + logsum(c); // + logsum(c); + by thchiu
	    //else P.w = P.w + logsum(c); //+ by thchiu
	    //return (newP);
		if( depth > P.depth) P.depth +=1; 
	}
	

	
#ifdef PROFILE
    if(para.dim == 1)
        timer.TimerFinish_cpu(timer.tTime_rand_choose_C , timer.sTime_rand_choose_C);
    else
        timer.TimerFinish_cpu(timer.tTime_rand_choose_NC , timer.sTime_rand_choose_NC);
#endif
}


bool sample_trees(vector<OnePartition_data>& Ps, vector<vector<double> >& c_originals, vector<int>&cutregs,  double maxpercentage, double* time_array, bool discrete, parameters& para, double beta, int Cdepth, string method, double tree_temperature, bool UQ_flag) {   //+ by thchiu
    vector<double> c;
    double smooth = 0.0;
    //OnePartition_data newP;   //- by thchiu

    for (int i = 0; i < (int) Ps.size(); i++) {
        //get_info_from_onetree(Ps[i], c_originals[i], c, cutregs[i], maxpercentage, smooth, time_array, discrete); //- by thchiu
        get_info_from_onetree(Ps[i], c_originals[i], c, cutregs[i], maxpercentage, smooth, time_array, discrete, para, Cdepth, method);   //+ by thchiu

        if (max(c) < 1 - MAXNUM) {
            return true;
        }

        //rand_choose_one_nextleveltree(newP, Ps[i], c,cutregs[i], smooth, discrete); //- by thchiu
        rand_choose_one_nextleveltree(Ps[i], c_originals[i], c,cutregs[i], smooth, discrete, para, beta, Cdepth, method,tree_temperature, UQ_flag); //+ by thchiu

        //Ps[i].sregs = newP.sregs;   //problem6: why? data copy, - by thchiu
        //Ps[i].w = newP.w; //- by thchiu
    }
	// get rid of normalization, as we hope weights are comparable across different trees under different perturbations, linxi, 9/10/24
    double largestwt = Ps[0].w;
    for(int i=0; i<(int)Ps.size(); i++)  {  //problem5: start from 1 is OK, thchiu
        if(Ps[i].w>largestwt)  largestwt= Ps[i].w;
    }
    for (int i = 0; i < (int) Ps.size(); i++) { // normalization, avoid too small or large weight
        Ps[i].w = Ps[i].w - largestwt + 5; 

    }
	
    return false;
}



// get the best partition so far: the one with the largest weight*number of appearances in the samplings.

// return the maximum of the score among the partitions in the level; modelBP stores the index of the maximum
double getmaxlBP(vector<OnePartition_data>& newPs, int & modelBP) {

    double maxlBP = newPs[0].post_p; //lBPscore(newPs[0]);
    modelBP=0;
    for (int i = 1; i < (int) newPs.size(); i++) {
		// rewrite for OPT, as we do not calculate posterior probability separately, linxi
		if(newPs[i].post_p > maxlBP) {
			maxlBP= newPs[i].post_p;
			modelBP=i;
		}
		/*
            if(lBPscore(newPs[i])>maxlBP){
                 maxlBP = lBPscore(newPs[i]);
                 modelBP=i;
            }
		*/
    }
    return (maxlBP);
}

// get the best partition so far: the one with the largest weight*number of appearances in the samplings.
// used for testing only, no need to run in the real data until the last run
void getmax_corr(vector<OnePartition_data>& newPs,  vector<int> & modeinds, vector<double>& modewts, double & ess) {

    vector<double> wts;
    //int dim = newPs[0].sregs[0].reg_code.size();  //commented by thchiu
    for (int i = 0; i < (int) newPs.size(); i++) { // normalization, avoid too small or large weight
        //       newPs[i].w = newPs[i].w - largestwt + 5;
        if(isnan(exp(newPs[i].w)))     wts.push_back(0.000001);
        else  wts.push_back(exp(newPs[i].w));

    }
    CompairSRegs_data mcompare0;
    sort(newPs.begin(), newPs.end(), mcompare0 );
    CompairSRegs_data mcompair;
    vector<double> wtprodtimes;
    int diffpartitionuntilnow=0;             // also the index of wtprodtimes;
    wtprodtimes.push_back( wts[0] ); //mouse
    modeinds.clear();
    modewts.clear();
    modeinds.push_back(0);

    for (int i = 1; i < (int) newPs.size(); i++) {
        if (!mcompair(newPs[i], newPs[i-1]) && !mcompair(newPs[i-1], newPs[i])) {
            wtprodtimes[diffpartitionuntilnow] += wts[i];
        } else {
            wtprodtimes.push_back(wts[i]);
            modeinds.push_back(i);
            diffpartitionuntilnow++;
        }

    }
    vector<pair<double, double> > prodweightsinfo;
    pair<double, double> onepair;
    for(int i=0; i<=diffpartitionuntilnow; i++){
        onepair.first = modeinds[i];
        onepair.second = wtprodtimes[i];
        prodweightsinfo.push_back(onepair);

    }
    sort(prodweightsinfo.begin(), prodweightsinfo.end(), myfunction4);
    modeinds.clear();
    modewts.clear();
    double largestweight = prodweightsinfo[0].second;
    bool more = true;
    for (int i = 0; more && i <= diffpartitionuntilnow; i++) {
        if ((prodweightsinfo[i].second / largestweight) > 0.0001) {
            modeinds.push_back((int)prodweightsinfo[i].first);
            modewts.push_back(prodweightsinfo[i].second);
        } else more = false;
    }

    cout<<"includeweight.size()="<<(int)modeinds.size()<<'\n';

    return;

}



//return the log of number of paths
/*
   double path_log_count(OnePartition & P) {
   cerr<<"count!\n";
   if (P.sregs.size() <= 3) return 0;
   int dim = P.sregs[0].reg_code.size();


   vector<int> avail_dim(dim, 1);
   for (int i = 0; i < dim; i++) {
   for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
   if (P.sregs[bri].reg_code[i].mask == 0) {
   avail_dim[i] = 0;
   break;
   }
   }
   }

   vector<double> log_count; // log(N) for each dimension, if this dimension is not avail, =0

   for (int i = 0; i < dim; i++) {
   if (avail_dim[i] == 1) {
   OnePartition P1, P2;
   int total1 = -1;
   int total2 = -1;
   for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
   shrink_reg_data_short temp(P.sregs[bri]);
   if ((~(temp.reg_code[i].mask >> 1)) & temp.reg_code[i].x) { // start from 0, represents the first half
   temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
   temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
   P1.sregs.push_back(temp);
   total1 += 1;
   } else { // start from 1, the other half
   temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
   temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
   P2.sregs.push_back(temp);
   total2 += 1;
   }
   }
   log_count.push_back(path_log_count(P1) + path_log_count(P2) + stirling(total1 + total2) - stirling(total1) - stirling(total2));

   }
   }
   return (logsum(log_count));
   }
 */



double path_log_count(OnePartition_data & P, vector< map< OnePartition, double, CompairSRegs> >& pathctmap) {
    int Pregsize = (int) P.sregs.size();
    if (Pregsize <= 3) return 0;
    if (Pregsize <= (int) pathctmap.size()) {
        OnePartition Pwodata;
        for (int k = 0; k < Pregsize; k++) {
            shrink_reg sr;
            sr.reg_code = P.sregs[k].reg_code;
            sr.num = 0;
            Pwodata.sregs.push_back(sr);
        }
        Pwodata.w = 0;

        sort(Pwodata.sregs.begin(), Pwodata.sregs.end(), myfunction3);
        map<OnePartition, double, CompairSRegs>::iterator it;
        it = pathctmap[Pregsize - 1].find(Pwodata);

        //        for( it2 = pathctmap[Pregsize - 1].begin(); it2!= pathctmap[Pregsize - 1].end(); it2++){
        //            CompairSRegs mouse_compare;
        //            cout<<mouse_compare(Pwodata, it2->first)<<" "<< mouse_compare(it2->first, Pwodata)<<"\n";

        //        }

        if (it == pathctmap[Pregsize - 1].end()) {
            cout << "not found in the map!\n";
            cout << "Pregsize-1=" << Pregsize - 1 << "\n";
            cout << "pathctmap[Pregsize-1].size=" << pathctmap[Pregsize - 1].begin()->first.sregs.size() << "\n";
            cout << "to find:\n";
            print_partition(Pwodata);
            cout<<"last one:\n";
            print_partition(pathctmap[Pregsize - 1].begin()->first);
        } else return it->second;
    }
    int dim = P.sregs[0].reg_code.size();

    vector<int> avail_dim(dim, 1);
    for (int i = 0; i < dim; i++) {
        for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
            if (P.sregs[bri].reg_code[i].mask == 0) {
                avail_dim[i] = 0;
                break;
            }
        }
    }

    vector<double> log_count; // log(N) for each dimension, if this dimension is not avail, =0

    for (int i = 0; i < dim; i++) {
        if (avail_dim[i] == 1) {
            OnePartition_data P1, P2;
            int total1 = -1;
            int total2 = -1;
            for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
                shrink_reg_data_short temp = P.sregs[bri];
                if ((~(temp.reg_code[i].mask >> 1)) & temp.reg_code[i].x) { // start from 0, represents the first half
                    temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
                    temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
                    P1.sregs.push_back(temp);
                    total1 += 1;
                }
                else { // start from 1, the other half
                    temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
                    temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
                    P2.sregs.push_back(temp);
                    total2 += 1;
                }
            }
            log_count.push_back(path_log_count(P1, pathctmap) + path_log_count(P2, pathctmap) + stirling(total1 + total2) - stirling(total1) - stirling(total2));
            //			cout<<"dim="<<dim<<" log_count="<<log_count[log_count.size()-1]<<endl;
        }
    }
    return (logsum(log_count));
}

/*
//up=mid+half
// low=mid-half

void generate_map(shrink_reg_data sr, vector<double> mid, vector<double> half, int dim, map<vector<usint_mask>, shrink_reg_data, CompairSReg>& newmap, vector<vector<double> >& data, int samplesize, int masklen) {

vector<int> subsubsize(2, 0);
vector<vector<int> > subsize(dim, subsubsize);

vector<vector<double> > onedata;
vector<vector<vector<double> > > d1(2, onedata);
vector< vector<vector<vector<double> > > > subdata(dim, d1);

for (int i = 0; i < samplesize; i++) {
for (int d = 0; d < dim; d++) {
if (data[i][d] < mid[d]) {
subsize[d][0] += 1;
subdata[d][0].push_back(data[i]);
} else {
subsize[d][1] += 1;
subdata[d][1].push_back(data[i]);
}
}
}

for (int d = 0; d < dim; d++) {
shrink_reg_data sr0 = sr;
shrink_reg_data sr1 = sr;
sr0.reg_code[d].x = sr.reg_code[d].x << 1;
sr1.reg_code[d].x = (sr.reg_code[d].x << 1) + 1;

sr0.reg_code[d].mask = (sr.reg_code[d].mask << 1) + 1;
sr1.reg_code[d].mask = sr0.reg_code[d].mask;

vector<double> mid0 = mid;
vector<double> mid1 = mid;
vector<double> half1 = half;

mid0[d] -= half[d]/2;
mid1[d] += half[d]/2;
half1[d]/=2;



if(subsize[d][0]>1 && newmap.find(sr0.reg_code)==newmap.end()){
sr0.num=subsize[d][0];
newmap[sr0.reg_code] = sr0;
if(half1[d]>pow(0.5,masklen) )   generate_map ( sr0, mid0, half1, dim, newmap, subdata[d][0], subsize[d][0], masklen);
}
if(subsize[d][1]>1 && newmap.find(sr1.reg_code)==newmap.end()){
sr1.num=subsize[d][1];
newmap[sr1.reg_code] = sr1;
if(half1[d]>pow(0.5,masklen) )  generate_map ( sr1, mid1, half1, dim, newmap, subdata[d][1], subsize[d][1], masklen);
}

}
return;
}


 */

//void resample(vector<OnePartition_data> &Ps, int halfweight, vector<vector<double> >& c_originals, vector<int> &cutregs) {    //commented by thchiu
void resample(vector<OnePartition_data>* &PsPtr, parameters& p, vector<vector<double> >* &c_originalsPtr, vector<int>* &cutregsPtr, bool UQ_flag) {     //added by thchiu
#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerStart_cpu(timer.sTime_resample_C);
    else
        timer.TimerStart_cpu(timer.sTime_resample_NC);
#endif
    //////////////  Ps is changed to *PsPtr in this function, thchiu ///////////////
    int halfweight = p.resampling;
    int size = (int) (*PsPtr).size();
    vector<double> ec(size, 0);
    for (int i = 0; i < size; i++) {
        if (halfweight == 2) {
            (*PsPtr)[i].w /= 2.0;
            ec[i] = exp((*PsPtr)[i].w);
        } else {
            ec[i] = exp((*PsPtr)[i].w);
            (*PsPtr)[i].w = 1;
        }
    }
    //	cout<<endl;
    //vector<OnePartition_data> newPs;  //commented by thchiu
    vector<OnePartition_data>* newPsPtr;  //added by thchiu
    newPsPtr = new vector<OnePartition_data>; //added by thchiu

    //vector<vector<double> > newc_originals; //Problem12: optimized, thchiu    //- by thchiu
    //vector<int> newcutregs;   //- by thchiu
    vector<vector<double> >* newc_originalsPtr; //+ by thchiu
    newc_originalsPtr = new vector<vector<double> >;// + by thchiu
    vector<int>* newcutregsPtr; //+ by thchiu
    newcutregsPtr = new vector<int>; //+ by thchiu

    int* Newdata1D_Index_train;
    Newdata1D_Index_train = new int[p.samplesize * p.n];
	int* Newdata1D_Index_vali;
	Newdata1D_Index_vali = new int[p.valisize * p.n];
	int* Newdata1D_Index_infer;
	Newdata1D_Index_infer = new int[p.infersize * p.n];

    for (int i = 0; i < size; i++) {
        int index = rand_int(ec);
        //newPs.push_back((*PsPtr)[index]); //problem3: data copy?, thchiu, commented by thchiu
        //newc_originals.push_back(c_originalsP[index]); //problem3: data_copy?, - by thchiu
        //newcutregs.push_back(cutregs[index]); //- by thchiu
        (*newPsPtr).push_back((*PsPtr)[index]); //+ by thchiu
        (*newc_originalsPtr).push_back((*c_originalsPtr)[index]); //+ by thchiu
        (*newcutregsPtr).push_back((*cutregsPtr)[index]); //+ by thchiu
        (*newPsPtr)[i].PartitionID=i; //+ by thchiu

        for(int j=0; j<p.samplesize; j++){ //+ by thchiu
            Newdata1D_Index_train[p.samplesize * i + j] = p.pt_start[p.samplesize*index+j];
        }
		
		if(p.valisize >0) {
	        for(int j=0; j<p.valisize; j++){ //+ by thchiu
	            Newdata1D_Index_vali[p.valisize * i + j] = p.pt_start_vali[p.valisize*index+j];
	        }
		}
       	
		if(p.dim ==1 && UQ_flag) {
	        for(int j=0; j<p.infersize; j++){ //+ by thchiu
	            Newdata1D_Index_infer[p.infersize * i + j] = p.pt_start_infer[p.infersize*index+j];
	        }	
		}

        for(unsigned j=0; j<(*newPsPtr)[i].sregs.size(); j++){
            (*newPsPtr)[i].sregs[j].region_start = ((*PsPtr)[index].sregs[j].region_start%p.samplesize) + i*p.samplesize;
			if(p.valisize >0) {
				(*newPsPtr)[i].sregs[j].region_start_vali = ((*PsPtr)[index].sregs[j].region_start_vali%p.valisize) + i*p.valisize;			
			}
			if(p.dim ==1 && UQ_flag) {
				(*newPsPtr)[i].sregs[j].region_start_infer = ((*PsPtr)[index].sregs[j].region_start_infer%p.infersize) + i*p.infersize;
			}
        }
    }
    //Ps=newPs; // - by thchiu
    //c_originals=newc_originals;   //- by thchiu
    //cutregs=newcutregs;   //- by thchiu
    delete PsPtr;   //+ by ed520
    PsPtr = newPsPtr; // + by thchiu
    //delete c_originalsPtr;     //+ by ed520
    delete c_originalsPtr;  //+ by ed520
    c_originalsPtr = newc_originalsPtr; //+ by thchiu
    //delete cutregsPtr;      //+ by ed520
    delete cutregsPtr;  //+ by ed520
    cutregsPtr = newcutregsPtr; //+ by thchiu
    delete [] p.pt_start;   //+ by thchiu
    p.pt_start = Newdata1D_Index_train; //+ by thchiu
	if(p.valisize >0) {
		delete [] p.pt_start_vali;
		p.pt_start_vali = Newdata1D_Index_vali; //+ by thchiu	
	}	
	if(p.dim ==1 && UQ_flag) {
		delete [] p.pt_start_infer;
		p.pt_start_infer = Newdata1D_Index_infer; 
	}

#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerFinish_cpu(timer.tTime_resample_C , timer.sTime_resample_C);
    else
        timer.TimerFinish_cpu(timer.tTime_resample_NC , timer.sTime_resample_NC);
#endif
    return;
}

// n is the number of points in P
double densities_from_partition(vector<double> &testdata,  OnePartition_data &P, int n){

    int regnum = (int)P.sregs.size();
    int dim = (int) P.sregs[0].reg_code.size();


    pair<double, double> onerange;

    for (int i = 0; i < regnum; i++) {
        bool ind = true;
        for (int d = 0; d < dim; d++) {

            onerange = convert_ranges(P.sregs[i].reg_code[d]);
            if (testdata[d] > onerange.second || testdata[d] <= onerange.first) {
                ind = false;
                break;
            }
        }
        if (ind == true) {
            return ((double) P.sregs[i].num / n / exp(lprod_usint_mask(P.sregs[i].reg_code)));
        }
    }
    cout << "wrong densities_from_partition!" << endl;
    for (int d = 0; d < dim; d++) {
        cout<<testdata[d]<<" ";
    }
    exit(10);
    return -1;


}

// Get top K smallest edge length along each dimension
// if number of subregions is samller than K, then it's simply # of subregions
void get_smallest_edge_length(OnePartition_data &P, int K) {
    int regnum = (int)P.sregs.size();
    int dim = (int) P.sregs[0].reg_code.size();
	
	if(K > regnum) K = regnum;
	
	vector<vector<double> > edge_length;
	edge_length.resize(dim);
	
	for(int d=0; d< dim; d++) {
		edge_length[d].resize(K);
		for(int l =0; l< K; l++) {
			edge_length[d][l]=1.0;
		}	
		for(int i=0; i < regnum; i++) {
			double temp_length = 1/(double)(P.sregs[i].reg_code[d].mask+1);
			for(int l =0; l <K; l++) {
				if(temp_length < edge_length[d][l]) {
					//cout<<temp_length<<endl;
					if(l < K-1) {
						for(int m = K-1; m >l; m--) {
							edge_length[d][m] = edge_length[d][m-1];
						}
					}
					edge_length[d][l] = temp_length;
					break;
				}	
			}
		}
	}
	
	for(int d=0; d< dim; d++) {
		cout<<"smallest edge_length along dimension "<<d<<" is:"<<endl;
		for(int l=0; l <K; l++) {
			cout<<edge_length[d][l]<<"\t";
		}
		cout<<endl;
	}
	
}

//the density of a point from a mixture partition density.
// Ps.size()==wts.size(); wts need not to be normalized.

double densities_from_mixpartition(vector<double> &testdata, vector<OnePartition_data> &Ps, vector<double> wts, int n) {

    double den = 0;
    double sumwts = sum(wts);
    for (int i = 0; i < (int) Ps.size(); i++) den = den + densities_from_partition(testdata, Ps[i], n) * wts[i] / sumwts;
    return den;
}

double densities_from_mixpartition(vector<double> &testdata, vector<OnePartition_data> &Ps, vector<double> wts, int n, double lattic) {

    double den = 0;
    double sumwts = sum(wts);
    for (int i = 0; i < (int) Ps.size(); i++) den = den + one_density(Ps[i], testdata, lattic) * wts[i] / sumwts;
    return den;
}


//double densities_from_partition(vector<double> &testdata, vector<vector<double> > P) {  //- by ed520
double densities_from_partition(vector<double> &testdata, vector<vector<double> >& P) { //+ by ed520
    double map_den = -1.0;
    int dim = (int)testdata.size();
    bool found = false;
    //int visit = 0;	//+ 412 for testing
    for (int j = 0; !found && j < (int) P.size(); j++) {
        //visit++;
        found = true;
        for (int d = 0; found && d < dim; d++) {

            double low = P[j][2 * d];
            double high = P[j][(2 * d) + 1];

            if (testdata[d] <= low || testdata[d] > high) {
                found = false;
            }
        }

        if (found) {
            map_den = P[j][2 * dim];
        }
    }

    if (!found) {
        cerr << "Sample not found, idx:\n";
        map_den = 0;	//+ 412
    }
    //cout<<" Visit node amount: "<< visit << endl;	//+ 412 for testing
    //cout<< " PDF = "<< map_den << endl;
    return map_den;
}

double densities_from_tree(vector<double> &testdata, vector<vector<double> >& P){//+ 412
    double 	map_den = -1.0;
    int 	dim = (int)testdata.size();
    int 	Cur_ID = 0;
    bool 	found = false;
    bool 	leaf  = false;
    
    //cout << testdata[0] <<" "<< testdata[1] <<endl;
    
    while ( Cur_ID < (int) P.size() ){
        //cout<< Cur_ID <<endl;//
        leaf = false;
        if (P[Cur_ID][2*dim + 3] == -1){// If cut_dim = -1, it's a leaf node
            leaf = true;
            //cout<<" Leaf "<<endl;
        }        
        found = true;
        for (int d = 0; found && d < dim; d++){
            if (testdata[d] <= P[Cur_ID][2*d + 1] || testdata[d] > P[Cur_ID][2*d + 2]){
                found = false;
                //cout<<" Miss @ Node "<< Cur_ID <<": ";
            }
        }
        if (found && leaf){
            map_den = P[Cur_ID][2*dim + 1];
            //cout<<" Hit! "<<endl;
            break;
        }
        if (!found && (Cur_ID == 0)){// Data point is outside the root region
            cerr<<"Sample not found @ Root: "<<endl;
            map_den = 0;
            break;
        }        
        if (!found){// If the data point is not in this node, then stop and jump to other branch
            Cur_ID = DFS_NextNode_PDF(Cur_ID, P, dim);
        }
        else{
            Cur_ID++;
        }
    }
    //cout<<" PDF = "<<map_den<<endl<<endl;
    return map_den; 
}

int DFS_NextNode_PDF(int Cur_ID, vector<vector<double> >& P, int dim){ //+ 412; change name?
    int Next_ID = Cur_ID;
    int Parent_ID = -1;
    bool find_next = false;
    //cout << "Enter DFS_NxtNode" << endl;    
    while(!find_next){
        //cout<<" Cur_ID = "<< Next_ID <<" "<<" Par_ID = "<< P[Next_ID][2*dim + 4] <<endl;
        Parent_ID = P[Next_ID][2*dim + 4];
        if (Parent_ID == -1){
            break;
        }
        int Bro_ID = (P[Parent_ID][2*dim + 5] == Next_ID) ? P[Parent_ID][2*dim + 6] : P[Parent_ID][2*dim + 5];
        if (Bro_ID > Cur_ID){
            Next_ID = Bro_ID;
            find_next = true;
        }
        else{
            Next_ID = Parent_ID;
        }
        if (Next_ID == 0){// modified
            break;
        }
    }
    return Next_ID;
}

double CDF_from_partition(vector<double> &testdata, vector<vector<double> >& P) { //+ 412
    double map_den = -1.0;
    double prob_reg;			
    double length;			
    double volume;			
    double CDF = 0;			
    int    dim = (int)testdata.size();  
    bool   count_flag = false;			//+ 412; alternative mode for accumulating CDF, can be trun-off
    int    total_amount = 0;
    for (int j = 0; j < (int) P.size(); j++){
        total_amount += P[j][2*dim+1];
    }
    
    for (int j = 0; j < (int) P.size(); j++) { 	// for each region
        count_flag = true;
        prob_reg = 0;
        volume = 1;
        for (int d = 0; d < dim; d++) {		// for each dimension
            double low = P[j][2 * d];		// lower bound of j-th region @ the d-th dimension
            double high = P[j][2 * d + 1];	// upper bound ...       
            length = testdata[d] - low;
            if (length < 0){
                length = 0;
            }
            if (length > (high - low)){
                length = high - low;
            }
            volume *= length;
            count_flag = count_flag && (testdata[d] >= high);
        }
        if (count_flag){
            CDF += P[j][2 * dim + 1]/total_amount;
        }
        else{
            map_den = P[j][2 * dim];
            prob_reg = map_den*volume;
            CDF += prob_reg;
        }
    }
    return CDF;
}

double CDF_from_tree(vector<double> &testdata, vector<vector<double> >& P) { //+ 412
    double 	CDF = 0;
    double	length;
    int 	dim = (int)testdata.size();
    int 	Cur_ID = 0;
    bool 	leaf  = false;
    bool	cover = false;
    //cout << "Enter CDF_from_tree." <<endl;    
    while ( Cur_ID < (int) P.size() ){
        //cout<< Cur_ID <<"->";
        leaf = false;
        if (P[Cur_ID][2*dim + 3] == -1){// If cut_dim = -1, it's a leaf node
            leaf = true;
            //cout<<" Leaf "<<endl;
        }        
        cover = true;
        double	volume = 1;
        for (int d = 0; d < dim; d++){
            double low = P[Cur_ID][2*d + 1];
            double high = P[Cur_ID][2*d + 2];
            length = testdata[d] - low;
            if (testdata[d] <= low){
                length = 0;
            }
            if (testdata[d] > high){
                length = high - low;                
            }
            cover = cover && (testdata[d] >= high);
            volume *= length;
        }
        //cout<<" Cover: "<< cover <<" Leaf: "<< leaf <<" Vol = "<< volume <<" ;  ";           
        if (cover){
            CDF += P[Cur_ID][2*dim + 2]/P[0][2*dim + 2];
            Cur_ID = DFS_NextNode_PDF(Cur_ID, P, dim); // Since all children cannot be visited, finding another branch
            //cout<<"(I)";
            if (Cur_ID == 0){
                break;
            }
        }
        else if ( leaf && (volume > 0) ){
            CDF += P[Cur_ID][2*dim + 1] * volume;
            Cur_ID++;
            //cout<<"II"<<endl;
        }
        else if ( volume == 0 ){
            Cur_ID = DFS_NextNode_PDF(Cur_ID, P, dim); // Since all children cannot be visited, finding another branch
            //cout<<"(III)";
            if (Cur_ID == 0){
                break;
            }
        }
        else{
            Cur_ID++;
            //cout<<"(IV)";
        }
        //cout<<Cur_ID<<"->";
    }
    //cout<<endl;
    return CDF; 
}	
