#include "density_tree_forest.h"

void bma_density_forest(parameters & p, parameters & para, vector<OneDensity_data>& ensemble, int ensemble_sizeP, int n_shiftM, double shift_sizeM, int n_shiftP, double shift_sizeP, double betaM, double betaP, string ofilename) {
	int dim = p.dim;
	
	// Release the memory 
	delete[] p.data1D;
	delete[] p.data1D_vali;
	//delete[] p.data1D_infer;
  	delete[] p.pt_start;
	delete[] p.pt_start_vali;
	//delete[] p.pt_start_infer;
		
	p.data1D = nullptr;
	p.data1D_vali = nullptr;
	//p.data1D_infer = nullptr;
	p.pt_start = nullptr;
	p.pt_start_vali = nullptr;
	//p.pt_start_infer = nullptr;

    p.data1D= new double[p.dim*p.samplesize];  // new
	p.data1D_vali = new double[p.dim*p.valisize]; 
    p.pt_start = new int [p.samplesize * p.n];
	p.pt_start_vali = new int [p.valisize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.valisize; j++){
            p.pt_start_vali[i*p.valisize+j] = j*p.dim;
        }
    }
	
	delete[] para.data1D;
	delete[] para.data1D_vali;
	//delete[] para.data1D_infer;
	delete[] para.pt_start;
	delete[] para.pt_start_vali;
	//delete[] para.pt_start_infer;

	para.data1D = nullptr;
	para.data1D_vali = nullptr;
	//para.data1D_infer = nullptr;
	para.pt_start = nullptr;
	para.pt_start_vali = nullptr;
	//para.pt_start_infer = nullptr;

    para.data1D = new double[para.dim*para.samplesize]; // new
	para.data1D_vali = new double[para.dim*para.valisize];
	if(p.dim == 1) {
		//para.data1D_infer = new double[para.dim*para.infersize];
	}
    para.pt_start = new int[para.samplesize*para.n];
	para.pt_start_vali = new int[para.valisize*para.n];
	if(p.dim ==1) {
		//para.pt_start_infer = new int[para.infersize*para.n];
	}
    //    bool pathct=false;
    //    bool discrete=true;

    //    double smoothneighbordist=0.001;
    // density to be learned at the current step, linxi 
	//   vector<OnePartition_data>  bestlevelPs;
    OnePartition_data  bestlevelPs;
    //    int bestind;
	
	ensemble.clear();
	
	unsigned seed = time(0);
	srand(seed);
	
	for(int l=0; l< ensemble_sizeP; l++) {
		cout<<"-----------------------------"<<endl;
		cout<<"------------l="<<l<<"--------------"<<endl;
		cout<<"-----------------------------"<<endl;
		// Output related, deleted to make the computation time shorter.
		/*
		string tempdir = ofilename+"BDF/shift_"+std::to_string(l)+"/";
		string cline = "mkdir -p "+tempdir;
		const int dir= system(cline.c_str());
		if (dir< 0) return;
		*/
		OneDensity_data onefit;
		onefit.shift.clear();
		onefit.joint_shift.clear();
		vector<int> shift_int(p.dim, 0);
		
		// Output related, deleted to make the computation time shorter.
		//string tempshift = tempdir+"Mshift.txt";
		//ofstream outfileMshift(tempshift.c_str());
		for(int d=0; d <p.dim; d++) {
			shift_int[d] = rand_int(0, (n_shiftM-1));
			onefit.shift.push_back( (double)shift_int[d] * shift_sizeM );
			cout<<"For marginal distribution, shift along dimension "<<d<<" is: "<<onefit.shift[d]<<"\t"<<shift_int[d]<<endl;
			//outfileMshift<<onefit.shift[d]<<endl;
		}
		//outfileMshift.close();
	
		//tempshift = tempdir +"joint_shift.txt";
		//ofstream outfilePshift(tempshift.c_str());
		for(int d=0; d<p.dim; d++) {
			shift_int[d] = rand_int(0, (n_shiftP-1));
			onefit.joint_shift.push_back( (double)shift_int[d] * shift_sizeP);
			cout<<"For joint distribution, shift along dimension "<<d<<" is: "<<onefit.joint_shift[d]<<"\t"<<shift_int[d]<<endl;
			//outfilePshift<<onefit.joint_shift[d]<<endl;
		}
		//outfilePshift.close();
	
		// Output related, deleted to make the computation time shorter.
		//string temptrans = tempdir+"trans_data.txt"; //file for Copula transformed data
		//cout<<temptrans<<endl;
		
		onefit.marginalP.clear();
	    //ofstream outfiletrans (temptrans.c_str());
	    for(int d=0; d<dim; d++){
	        cerr<<"dim="<<d<<endl<<endl;
	        para.data.clear();
			para.validata.clear();
	        vector<double> onepiece(1,0); //initialization for a vector of doubles, linxi
	        vector<double> onevariable(para.samplesize, 0);
			vector<double> onevariable_vali(para.valisize, 0);
			//vector<double> oneshift(1,0);
	
	
	        // shift the data along the dth dimension, linxi
			// first shift the original data, 
			// then each coordinate is divided by 2, linxi
			// update both para.data and para.data1D, linxi
			// p.data remains unchanged, linxi
			for(int i=0; i<para.samplesize; i++)   {
	            onepiece[0]=(p.data[i][d]+onefit.shift[d])*0.5;
	            onevariable[i]= (p.data[i][d]+onefit.shift[d])*0.5;
	            para.data.push_back(onepiece);

	            para.data1D[i*para.dim]=(p.data[i][d]+onefit.shift[d])*0.5; // new
	            for(int j=0; j<para.n; j++)
	                para.pt_start[j*para.samplesize+i] = i*para.dim;
	        }
			
	        for(int i=0; i<para.valisize; i++)   {
	            onepiece[0]=(p.validata[i][d]+onefit.shift[d])*0.5;
	            onevariable_vali[i]= (p.validata[i][d]+onefit.shift[d])*0.5;
	            para.validata.push_back(onepiece);

	            para.data1D_vali[i*para.dim]=(p.validata[i][d]+onefit.shift[d])*0.5; // new
	            for(int j=0; j<para.n; j++)
	                para.pt_start_vali[j*para.valisize+i] = i*para.dim;
	        }
	

			// fit BSP for marginal dist, linxi
			// for bma, each tree is grown under the opt prior, without data splitting
			// "30, 100" are large tree depth limit 
	        SIS(para, true, bestlevelPs, betaM, 30, 100, "opt");
	
			// print top 5 smallest edge length;
			if(l/50==0) get_smallest_edge_length(bestlevelPs, 5);

	        onefit.marginalP.push_back( bestlevelPs);
			//    print_partition(bestlevelPs);
	        print_partition(onefit.marginalP[d]); // output the marginal dist to terminal, linxi

	        // get marginal quantiles (copula), linxi
	        vector<double> newx = Ftransform(onevariable, onefit.marginalP[d]); // defined in output.h
			// still scale to [0, 0.5] after marginal copula
			scale_by_two(&newx[0], para.samplesize, para.dim, true);
	
			// checkings, linxi
			cerr<<"orginal data is: ";
	        for(int i=0; i<10; i++) cerr<<p.data[i][d]<<" ";
	        cerr<<endl;

	        cerr<<"transformed data is: ";
			for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
	        cerr<<endl;

			// checkings, inverse transformation, linxi
			scale_by_two(&newx[0], para.samplesize, para.dim, false);
	        vector<double> recoveronevariable = inv_Ftransform(newx, onefit.marginalP[d]);
			scale_by_two(&recoveronevariable[0], para.samplesize, para.dim, false);

			cerr<<"recoverd data is: ";
	        for(int i=0; i<10; i++) cerr<<recoveronevariable[i]-onefit.shift[d]<<" ";
	        cerr<<endl;
	
			// update parameter set p for multivariate density estimation
			// updata the data component with the one after copula transformation, linxi
			// output the transformed data
			// transformed data is in [0,1], linxi
	        for(int i=0; i<p.samplesize; i++){
	            //dataF[i][d]=newx[i];            //new
	            //p.data[i][d] = newx[i];       //- by ed520, we use below
	            p.data1D[i*p.dim+d] = newx[i];  //+ by ed520
				// Output related, deleted to make the computation time shorter.
				//if(i< 100) outfiletrans << newx[i]<<'\t';
	        }
			
			vector<double> newx_vali = Ftransform(onevariable_vali, onefit.marginalP[d]);
		
	        for(int i=0; i<p.valisize; i++){
	            //dataF[i][d]=newx[i];            //new
	            //p.data[i][d] = newx[i];       //- by ed520, we use below
	            p.data1D_vali[i*p.dim+d] = newx_vali[i];  //+ by ed520
				// Output related, deleted to make the computation time shorter.
	            //if(i<100) outfiletrans << newx_vali[i]<<'\t';
	        }
			// Output related, deleted to make the computation time shorter.
	        //outfiletrans<<'\n';
		}
		// Output related, deleted to make the computation time shorter.
		//outfiletrans.close();
	
		if(p.dim > 1) {
			//cout<<"onefit.joint_shift[0]: "<<onefit.joint_shift[0]<<" "<<onefit.joint_shift[1]<<endl;
			//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
			shift_data(p.data1D, p.samplesize, p.dim, onefit.joint_shift);
			//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
			scale_by_two(p.data1D, p.samplesize, p.dim, true);
			
			shift_data(p.data1D_vali, p.valisize, p.dim, onefit.joint_shift);
			scale_by_two(p.data1D_vali, p.valisize, p.dim, true);
			//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
		    //    para.dim = dim;
		    //    para.data.clear();
		    //    para.data=dataF;
		    //    para.levels = p.levels;
		    //    para.n = p.n;
		    //    para.steps = p.steps;
		    //    para.resampling = p.resampling;          //new
		    //   SIS(para, true, bestlevelPs, bestind);

		    // fit the multivariate BSP after copula trans, linxi
			//onefit.jointP.clear();
			SIS(p, true, onefit.jointP, betaP, 30, 100, "opt");

			// print top 5 smallest edge length;
			if(l/100 ==0) get_smallest_edge_length(onefit.jointP, 5);

			// output marginal dist's to files, linxi
		    //     print_Finvtransform(P, marginalP);
			// Output related, deleted to make the computation time shorter.
			/*
		    for (int d = 0; d < dim; d++) {
		        string temp= tempdir+"dim"+std::to_string(d)+".txt";
		        ofstream outfile (temp.c_str());

		        // fprintf(pFile, "d= %d\n", d);
		        print_partition(outfile, onefit.marginalP[d]);
		        outfile.close();
		    }
			*/
			
			// Output related, deleted to make the computation time shorter.
			/*
			string distTrans= tempdir+"transJoint.txt";
			string distOrig=tempdir+"origJoint.txt";	//+ 412
			ofstream outdistTrans(distTrans.c_str());
			ofstream outdistOrig(distOrig.c_str());		//+ 412
			*/
			OnePartition_data &op = onefit.jointP;		

			//    outfile <<endl << "After transform back" << endl;
			for (int i = 0; i < (int) op.sregs.size(); i++) {
				pair<double, double> range;
			    for (int d = 0; d < dim; d++) {
					range = convert_ranges(op.sregs[i].reg_code[d]);
			        double ub,lb;
			        ub = inv_Ftransform(range.second*2-onefit.joint_shift[d], onefit.marginalP[d]);
			        // ub = ub * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
					ub = ub * 2.0 -onefit.shift[d] - 0.000001;
			        lb =  inv_Ftransform(range.first*2-onefit.joint_shift[d], onefit.marginalP[d]);
			        // lb = lb * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
					lb = lb *2.0 - onefit.shift[d] - 0.000001;
			        //         outfile<<range.first<<" "<<range.second<<" ";
					// Output related, deleted to make the computation time shorter.
			        //outdistOrig <<lb << " " << ub <<" ";
			        //outdistTrans << range.first << " " << range.second << " ";
				}
				// calculate the density on a subregion, linxi
				// Output related, deleted to make the computation time shorter.
			    //outdistTrans<< (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;
			  	//outdistOrig << (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)+(double)p.dim*log(2.0)) <<" "<< op.sregs[i].num<<endl;
			}
			
			// Output related, deleted to make the computation time shorter.
			//outdistTrans.close();
			//outdistOrig.close();	
		}else {
			// Output related, deleted to make the computation time shorter.
			/*
			string distOrig=tempdir+"origJoint.txt";	//+ 412
			ofstream outdistOrig(distOrig.c_str());	
			print_partition(outdistOrig, onefit.marginalP[0]);
			outdistOrig.close();	
			*/
		}
		// output joint dist to files, linxi	
	
		ensemble.push_back(onefit);		
	}
}

// int n_shiftM, double shift_sizeM, int n_shiftP, double shift_sizeP,
void gibbs_density_forest(parameters & p, parameters & para, vector<OneForest_data>& forest_of_Fs, int depth_max, int depth_min, int ensemble_sizeP, int number_of_forests, double betaM, double betaP, double splitP, string ofilename, double tree_temperature, double forest_temperature, int Cdepth, bool UQ_flag) {
	int dim = p.dim;
	
	// Release the memory 
	delete[] p.data1D;
	delete[] p.data1D_vali;
	//delete[] p.data1D_infer;
  	delete[] p.pt_start;
	delete[] p.pt_start_vali;
	//delete[] p.pt_start_infer;
		
	p.data1D = nullptr;
	p.data1D_vali = nullptr;
	//p.data1D_infer = nullptr;
	p.pt_start = nullptr;
	p.pt_start_vali = nullptr;
	//p.pt_start_infer = nullptr;

    p.data1D= new double[p.dim*p.samplesize];  // new
	p.data1D_vali = new double[p.dim*p.valisize]; 
    p.pt_start = new int [p.samplesize * p.n];
	p.pt_start_vali = new int [p.valisize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.valisize; j++){
            p.pt_start_vali[i*p.valisize+j] = j*p.dim;
        }
    }
	
	delete[] para.data1D;
	delete[] para.data1D_vali;
	if(p.dim ==1 && UQ_flag) {
		delete[] para.data1D_infer;
	}
	delete[] para.pt_start;
	delete[] para.pt_start_vali;
	if(p.dim ==1 && UQ_flag) {
		delete[] para.pt_start_infer;	
	}

	para.data1D = nullptr;
	para.data1D_vali = nullptr;
	if(p.dim == 1 && UQ_flag) {
		para.data1D_infer = nullptr;	
	}
	para.pt_start = nullptr;
	para.pt_start_vali = nullptr;
	if(p.dim ==1 && UQ_flag) {
		para.pt_start_infer = nullptr;	
	}

    para.data1D = new double[para.dim*para.samplesize]; // new
	para.data1D_vali = new double[para.dim*para.valisize];
	if(p.dim == 1 && UQ_flag) {
		para.data1D_infer = new double[para.dim*para.infersize];
	}
    para.pt_start = new int[para.samplesize*para.n];
	para.pt_start_vali = new int[para.valisize*para.n];
	if(p.dim ==1 && UQ_flag) {
		para.pt_start_infer = new int[para.infersize*para.n];
	}
	

    //    bool pathct=false;
    //    bool discrete=true;

    //    double smoothneighbordist=0.001;
    // density to be learned at the current step, linxi 
	//   vector<OnePartition_data>  bestlevelPs;
    OnePartition_data  bestlevelPs;
	
	//for(int eta =-1; eta < 2; eta++) {
	for(int b=0; b< number_of_forests; b++ ) {
		OneForest_data ensemble;
		ensemble.forestP.clear();
		ensemble.weight =0.0;
		
		// randomly sample the tree depth from the uniform dist on interval [depth_min, depth_max]
		int k = rand_int(depth_min, depth_max);
		
		//double tempNum = exp((double)(k+eta)*log(2.0));
		double tempNum = exp((double)(k)*log(2.0));
		int n_shift = (int)tempNum;
		//cout<<tempNum<<"\t"<<(double)n_shift +1 - tempNum <<"\t"<<tempNum - (double)n_shift<<endl;
		if((double)n_shift +1 - tempNum < tempNum - (double)n_shift) n_shift +=1;
			
		//int ensemble_size = ensemble_sizeP; // fix ensemble size
		//int ensemble_size = n_shift *p.dim*2;
		//if(ensemble_size > ensemble_sizeP) ensemble_size = ensemble_sizeP;
		//double grid_size = exp(-(double)(k+eta)*log(2.0));
		double grid_size = exp(-(double)(k)*log(2.0));
		double shift_size = grid_size * grid_size;
		for(int l=0; l< ensemble_sizeP; l++) {
			cout<<"-------------------------------------------"<<endl;
			cout<<"------------k="<<k<<"; l="<<l<<";--------------"<<endl;//"; eta="<<eta<<
			cout<<"-------------------------------------------"<<endl;
			/*
			string tempdir = ofilename+"GDF/forest_"+std::to_string(b)+"depth_"+std::to_string(k)+"/"+"shift_"+std::to_string(l)+"/"; 
			ofilename+"forest_"+std::to_string(b)+"depth_"+std::to_string(k)+"shift_"+std::to_string(l)+"/";
			string cline = "mkdir -p "+tempdir;
			const int dir= system(cline.c_str());
			if (dir< 0) return;
			*/
			OneDensity_data onefit;
			onefit.shift.clear();
			onefit.joint_shift.clear();
			vector<int> shift_int(p.dim, 0);

			
			//string tempshift = tempdir+"Mshift.txt";
			//ofstream outfileMshift(tempshift.c_str());
			for(int d=0; d <p.dim; d++) {
				shift_int[d] = rand_int(0, (n_shift-1));
				onefit.shift.push_back( (double)shift_int[d] * shift_size );
				cout<<"For marginal distribution, shift along dimension "<<d<<" is: "<<onefit.shift[d]<<"\t"<<shift_int[d]<<endl;
				//outfileMshift<<onefit.shift[d]<<endl;
			}
			//outfileMshift.close();
			

			
			//tempshift = tempdir +"joint_shift.txt";
			//ofstream outfilePshift(tempshift.c_str());
			// Do not apply the shift when we fit the joint dist
			for(int d=0; d<p.dim; d++) {
				shift_int[d] = 0;//rand_int(0, (n_shiftP-1));
				onefit.joint_shift.push_back( (double)shift_int[d] * shift_size);
				cout<<"For joint distribution, shift along dimension "<<d<<" is: "<<onefit.joint_shift[d]<<"\t"<<shift_int[d]<<endl;
				//outfilePshift<<onefit.joint_shift[d]<<endl;
			}
			//outfilePshift.close();
			

			// Output related, deleted.
			// string temptrans = tempdir+"trans_data.txt"; //file for Copula transformed data
			//cout<<temptrans<<endl;
	
			onefit.marginalP.clear();
			// Output related, deleted.
			// ofstream outfiletrans (temptrans.c_str());
			for(int d=0; d<p.dim; d++){
				cerr<<"dim="<<d<<endl<<endl;
				para.data.clear();
				para.validata.clear();
				if(UQ_flag) {
					para.inferdata.clear();
				}
				vector<double> onepiece(1,0); //initialization for a vector of doubles, linxi
				vector<double> onevariable(para.samplesize, 0);
				vector<double> onevariable_vali(para.valisize, 0);
				//vector<double> oneshift(1,0);


				// shift the data along the dth dimension, linxi
				// first shift the original data, 
				// then each coordinate is divided by 2, linxi
				// modify both para.data and para.data1D, linxi
				// p.data remains unchanged, linxi
				for(int i=0; i<para.samplesize; i++)   {
					onepiece[0]=(p.data[i][d]+onefit.shift[d])*0.5;
					onevariable[i]= (p.data[i][d]+onefit.shift[d])*0.5;
					para.data.push_back(onepiece);

					para.data1D[i*para.dim]=(p.data[i][d]+onefit.shift[d])*0.5; // new
					for(int j=0; j<para.n; j++)
						para.pt_start[j*para.samplesize+i] = i*para.dim;
				}
		
				for(int i=0; i<para.valisize; i++)   {
					onepiece[0]=(p.validata[i][d]+onefit.shift[d])*0.5;
					onevariable_vali[i]= (p.validata[i][d]+onefit.shift[d])*0.5;
					para.validata.push_back(onepiece);

					para.data1D_vali[i*para.dim]=(p.validata[i][d]+onefit.shift[d])*0.5; // new
					for(int j=0; j<para.n; j++)
						para.pt_start_vali[j*para.valisize+i] = i*para.dim;
				}
				
				if(p.dim ==1 && UQ_flag) {
					for(int i=0; i<para.infersize; i++)   {
						onepiece[0]=(p.inferdata[i+l*para.infersize][d]+onefit.shift[d])*0.5;
						para.inferdata.push_back(onepiece);

						para.data1D_infer[i*para.dim]=(p.inferdata[i+l*para.infersize][d]+onefit.shift[d])*0.5; // new
						for(int j=0; j<para.n; j++)
							para.pt_start_infer[j*para.infersize+i] = i*para.dim;
					}					
				}

				// fit BSP for marginal dist, linxi
				SIS(para, true, bestlevelPs, betaM, 30, 100, "gibbs", tree_temperature, UQ_flag);

				// print top 5 smallest edge length;
				if(l/50==0) get_smallest_edge_length(bestlevelPs, 5);

				onefit.marginalP.push_back( bestlevelPs);
				//    print_partition(bestlevelPs);
				print_partition(onefit.marginalP[d]); // output the marginal dist to terminal, linxi

				// get marginal quantiles (copula), linxi
				vector<double> newx = Ftransform(onevariable, onefit.marginalP[d]); // defined in output.h
				// still scale to [0, 0.5] after marginal copula
				scale_by_two(&newx[0], para.samplesize, para.dim, true);

				// checkings, linxi
				cerr<<"original data is: ";
				for(int i=0; i<10; i++) cerr<<p.data[i][d]<<" ";
				cerr<<endl;

				cerr<<"transformed data is: ";
				for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
				cerr<<endl;

				// checkings, inverse transformation, linxi
				scale_by_two(&newx[0], para.samplesize, para.dim, false);
				vector<double> recoveronevariable = inv_Ftransform(newx, onefit.marginalP[d]);
				scale_by_two(&recoveronevariable[0], para.samplesize, para.dim, false);

				cerr<<"recovered data is: ";
				for(int i=0; i<10; i++) cerr<<recoveronevariable[i]-onefit.shift[d]<<" ";
				cerr<<endl;

				// update parameter set p for multivariate density estimation
				// update the data component with the one after copula transformation, linxi
				// output the transformed data
				// transformed data is in [0,1], linxi
				// As we only perform UQ for 1-D examples, the following steps are not needed for inference data for multivariate distributions
				for(int i=0; i<p.samplesize; i++){
					//dataF[i][d]=newx[i];            //new
					//p.data[i][d] = newx[i];       //- by ed520, we use below
					p.data1D[i*p.dim+d] = newx[i];  //+ by ed520
					// Output related, deleted.
					// if(i< 100) outfiletrans << newx[i]<<'\t';
				}
		
				vector<double> newx_vali = Ftransform(onevariable_vali, onefit.marginalP[d]);
	
				for(int i=0; i<p.valisize; i++){
					//dataF[i][d]=newx[i];            //new
					//p.data[i][d] = newx[i];       //- by ed520, we use below
					p.data1D_vali[i*p.dim+d] = newx_vali[i];  //+ by ed520
					// Output related, deleted.
					// if(i<100) outfiletrans << newx_vali[i]<<'\t';
				}
				// Output related, deleted.
				// outfiletrans<<'\n';
			}

			// Output related, deleted.
			// outfiletrans.close();

			if(p.dim > 1) {
				//cout<<"onefit.joint_shift[0]: "<<onefit.joint_shift[0]<<" "<<onefit.joint_shift[1]<<endl;
				//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
				shift_data(p.data1D, p.samplesize, p.dim, onefit.joint_shift);
				//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
				scale_by_two(p.data1D, p.samplesize, p.dim, true);
		
				shift_data(p.data1D_vali, p.valisize, p.dim, onefit.joint_shift);
				scale_by_two(p.data1D_vali, p.valisize, p.dim, true);
				//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
				//    para.dim = dim;
				//    para.data.clear();
				//    para.data=dataF;
				//    para.levels = p.levels;
				//    para.n = p.n;
				//    para.steps = p.steps;
				//    para.resampling = p.resampling;          //new
				//   SIS(para, true, bestlevelPs, bestind);

				// fit the multivariate BSP after copula trans, linxi
				//onefit.jointP.clear();
				SIS(p, true, onefit.jointP, betaP, 30, 100, "gibbs", tree_temperature);

				// print top 5 smallest edge length;
				if(l/100 ==0) get_smallest_edge_length(onefit.jointP, 5);

				// output marginal dist's to files, linxi
				//     print_Finvtransform(P, marginalP);
				// Output related, deleted.
				/*
				for (int d = 0; d < dim; d++) {
					string temp= tempdir+"dim"+std::to_string(d)+".txt";
					ofstream outfile (temp.c_str());

					// fprintf(pFile, "d= %d\n", d);
					print_partition(outfile, onefit.marginalP[d]);
					outfile.close();
				}
				*/
	
				// Output related, deleted.
				/*
				string distTrans= tempdir+"transJoint.txt";
				string distOrig=tempdir+"origJoint.txt";	//+ 412
				ofstream outdistTrans(distTrans.c_str());
				ofstream outdistOrig(distOrig.c_str());		//+ 412
				*/
				OnePartition_data &op = onefit.jointP;		

				//    outfile <<endl << "After transform back" << endl;
				for (int i = 0; i <  (int) op.sregs.size(); i++) {
					pair<double, double> range;
					for (int d = 0; d < dim; d++) {
						range = convert_ranges(op.sregs[i].reg_code[d]);
						double ub,lb;
						ub = inv_Ftransform(range.second*2-onefit.joint_shift[d], onefit.marginalP[d]);
						// ub = ub * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
						ub = ub * 2.0 -onefit.shift[d] - 0.000001;
						lb =  inv_Ftransform(range.first*2-onefit.joint_shift[d], onefit.marginalP[d]);
						// lb = lb * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
						lb = lb *2.0 - onefit.shift[d] - 0.000001;
						//         outfile<<range.first<<" "<<range.second<<" ";
						// Output related, deleted.
						// outdistOrig <<lb << " " << ub <<" ";
						// Output related, deleted.
						// outdistTrans << range.first << " " << range.second << " ";
					}
					// calculate the density on a subregion, linxi
					// Output related, deleted.
					// outdistTrans<< (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;
					// outdistOrig << (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)+(double)p.dim*log(2.0)) <<" "<< op.sregs[i].num<<endl;
				}
				// Output related, deleted.
				// outdistTrans.close();
				// outdistOrig.close();	
			}else {
				// Output related, deleted.
				/*
				string distOrig=tempdir+"origJoint.txt";	//+ 412
				ofstream outdistOrig(distOrig.c_str());	
				print_partition(outdistOrig, onefit.marginalP[0]);
				outdistOrig.close();
				*/	
			}
			// output joint dist to files, linxi
		
			ensemble.forestP.push_back(onefit);		
		}
	
		vector<double> gibbs_weights = Gibbs_FOF_weight(p.validata, ensemble.forestP, betaM, betaP, splitP, tree_temperature, forest_temperature);
		ensemble.weight = gibbs_weights[0];
		ensemble.post_p = gibbs_weights[1];
		ensemble.tree_weight = gibbs_weights[2]; 
		//ensemble.weight -= 0.2* ((double)k )* exp((double)(k)*log(2.0));
		cout<<"The weight for forest "<<k<<" is:"<<ensemble.weight<<endl;
		forest_of_Fs.push_back(ensemble);					
	}
	//}
}

void optional_polya_tree_copula(parameters & p, parameters & para, vector<OnePartition_data>& marginalP, OnePartition_data & bestlevelPs, double betaM, double betaP, string ofilename){
	int dim = p.dim;
	
	// Release the memory 
	delete[] p.data1D;
	delete[] p.data1D_vali;
	//delete[] p.data1D_infer;
  	delete[] p.pt_start;
	delete[] p.pt_start_vali;
	//delete[] p.pt_start_infer;
		
	p.data1D = nullptr;
	p.data1D_vali = nullptr;
	//p.data1D_infer = nullptr;
	p.pt_start = nullptr;
	p.pt_start_vali = nullptr;
	//p.pt_start_infer = nullptr;

    p.data1D= new double[p.dim*p.samplesize];  // new
	p.data1D_vali = new double[p.dim*p.valisize]; 
    p.pt_start = new int [p.samplesize * p.n];
	p.pt_start_vali = new int [p.valisize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.valisize; j++){
            p.pt_start_vali[i*p.valisize+j] = j*p.dim;
        }
    }
	
	delete[] para.data1D;
	delete[] para.data1D_vali;
	//delete[] para.data1D_infer;
	delete[] para.pt_start;
	delete[] para.pt_start_vali;
	//delete[] para.pt_start_infer;

	para.data1D = nullptr;
	para.data1D_vali = nullptr;
	//para.data1D_infer = nullptr;
	para.pt_start = nullptr;
	para.pt_start_vali = nullptr;
	//para.pt_start_infer = nullptr;

    para.data1D = new double[para.dim*para.samplesize]; // new
	para.data1D_vali = new double[para.dim*para.valisize];
	if(p.dim == 1) {
		//para.data1D_infer = new double[para.dim*para.infersize];
	}
    para.pt_start = new int[para.samplesize*para.n];
	para.pt_start_vali = new int[para.valisize*para.n];
	if(p.dim ==1) {
		//para.pt_start_infer = new int[para.infersize*para.n];
	}
	
	// scale_by_two(p.data, false);
	string temptrans = ofilename+"OPT/trans_data.txt"; //file for Copula transformed data
    ofstream outfiletrans (temptrans.c_str());
    for(int d=0; d<dim; d++){
        cerr<<"dim="<<d<<endl<<endl;
        para.data.clear();
		para.validata.clear();
        vector<double> onepiece(1,0); //initialization for a vector of doubles, linxi
        vector<double> onevariable(para.samplesize, 0);
		vector<double> onevariable_vali(para.valisize, 0);
	
	
        for(int i=0; i<para.samplesize; i++)   {
            onepiece[0]=p.data[i][d];
            onevariable[i]= p.data[i][d];
            para.data.push_back(onepiece);

            para.data1D[i*para.dim]=p.data[i][d]; // new
            for(int j=0; j<para.n; j++)
                para.pt_start[j*para.samplesize+i] = i*para.dim;
        }
		
        for(int i=0; i<para.valisize; i++)   {
            onepiece[0]=p.validata[i][d];
            onevariable_vali[i]= p.validata[i][d];
            para.validata.push_back(onepiece);

            para.data1D_vali[i*para.dim]=p.validata[i][d]; // new
            for(int j=0; j<para.n; j++)
                para.pt_start_vali[j*para.valisize+i] = i*para.dim;
        }
	

		// fit BSP for marginal dist, linxi
        //SIS(para, true, bestlevelPs, betaM);
		SIS(para, true, bestlevelPs, betaM, 30, 100, "opt");
	
		// print top 5 smallest edge length;
		get_smallest_edge_length(bestlevelPs, 5);

        marginalP.push_back( bestlevelPs);
		//    print_partition(bestlevelPs);
        print_partition(marginalP[d]); // output the marginal dist to terminal, linxi

        // get marginal quantiles (copula), linxi
        vector<double> newx= Ftransform(onevariable, marginalP[d]); // defined in output.h
		//scale_by_two(&newx[0], para.samplesize, para.dim, true);
	
		// checkings, linxi
		cerr<<"orginal data is: ";
        for(int i=0; i<10; i++) cerr<<onevariable[i]<<" ";
        cerr<<endl;

        cerr<<"transformed data is: ";
		for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
        cerr<<endl;

		// checkings, inverse transformation, linxi
		//scale_by_two(&newx[0], para.samplesize, para.dim, false);
        vector<double> recoveronevariable = inv_Ftransform(newx, marginalP[d]);

		cerr<<"recoverd data is: ";
        for(int i=0; i<10; i++) cerr<<recoveronevariable[i]<<" ";
        cerr<<endl;
	
		// update parameter set p for multivariate density estimation
		// updata the data component with the one after copula transformation, linxi
        for(int i=0; i<para.samplesize; i++){
            //dataF[i][d]=newx[i];            //new
            //p.data[i][d] = newx[i];       //- by ed520, we use below
            p.data1D[i*p.dim+d] = newx[i];  //+ by ed520
            if(i < 100) outfiletrans << newx[i]<<'\t';
        }
		
		vector<double> newx_vali = Ftransform(onevariable_vali, marginalP[d]);
		
        for(int i=0; i<p.valisize; i++){
            //dataF[i][d]=newx[i];            //new
            //p.data[i][d] = newx[i];       //- by ed520, we use below
            p.data1D_vali[i*p.dim+d] = newx_vali[i];  //+ by ed520
            if(i < 100) outfiletrans << newx_vali[i]<<'\t';
        }
        outfiletrans<<'\n';
    }

    outfiletrans.close();
	
    // fit the multivariate BSP after copula trans, linxi
	//SIS(p, true, bestlevelPs, betaP);
	SIS(p, true, bestlevelPs, betaP, 30, 100, "opt");

	// print top 5 smallest edge length;
	get_smallest_edge_length(bestlevelPs, 5);

	// output marginal dist's to files, linxi
    //     print_Finvtransform(P, marginalP);
    for (int d = 0; d < dim; d++) {
        string temp= ofilename+"OPT/dim"+toStr<int>(d)+".txt";
        ofstream outfile (temp.c_str());

        // fprintf(pFile, "d= %d\n", d);
        print_partition(outfile, marginalP[d]);
        outfile.close();
    }

	// output joint dist to files, linxi
    string temp= ofilename+"OPT/transJoint.txt";
    string temp2=ofilename+"OPT/origJoint.txt";	//+ 412
    ofstream outfile (temp.c_str());
    ofstream outfile2(temp2.c_str());		//+ 412
    //   print_partition(outfile, bestlevelPs[bestind]);
    //   OnePartition_data op = bestlevelPs[bestind];
    OnePartition_data &op = bestlevelPs;

    //    outfile <<endl << "After transform back" << endl;
    for (int i = 0; i < (int) op.sregs.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < dim; d++) {
            range = convert_ranges(op.sregs[i].reg_code[d]);
            double ub,lb;
            ub = inv_Ftransform(range.second, marginalP[d]);
            //ub = ub * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
            lb =  inv_Ftransform(range.first, marginalP[d]);
            //lb = lb * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
            //         outfile<<range.first<<" "<<range.second<<" ";
            outfile2 <<lb << " " << ub <<" ";
            outfile << range.first << " " << range.second << " ";
        }
		// calculate the density on a subregion, linxi
        outfile<< (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;
        outfile2 << (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;
    }

    outfile.close();
    outfile2.close();	

    //+ by ed520
    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<p.dim<<endl;
    cout<<"Data_num = "<<p.samplesize+p.valisize<<endl;
    cout<<"# of Partition = "<<p.n<<endl;
}

void density_RF(parameters & p, parameters & para, vector<OneDensity_data>& ensemble, int ensemble_sizeP, int n_shiftM, double shift_sizeM, int n_shiftP, double shift_sizeP, double betaM, double betaP, string ofilename, int Cdepth) {
	int dim = p.dim;
	
	// Release the memory 
	delete[] p.data1D;
	delete[] p.data1D_vali;
	//delete[] p.data1D_infer;
  	delete[] p.pt_start;
	delete[] p.pt_start_vali;
	//delete[] p.pt_start_infer;
		
	p.data1D = nullptr;
	p.data1D_vali = nullptr;
	//p.data1D_infer = nullptr;
	p.pt_start = nullptr;
	p.pt_start_vali = nullptr;
	//p.pt_start_infer = nullptr;

    p.data1D= new double[p.dim*p.samplesize];  // new
	p.data1D_vali = new double[p.dim*p.valisize]; 
    p.pt_start = new int [p.samplesize * p.n];
	p.pt_start_vali = new int [p.valisize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.valisize; j++){
            p.pt_start_vali[i*p.valisize+j] = j*p.dim;
        }
    }
	
	delete[] para.data1D;
	delete[] para.data1D_vali;
	//delete[] para.data1D_infer;
	delete[] para.pt_start;
	delete[] para.pt_start_vali;
	//delete[] para.pt_start_infer;

	para.data1D = nullptr;
	para.data1D_vali = nullptr;
	//para.data1D_infer = nullptr;
	para.pt_start = nullptr;
	para.pt_start_vali = nullptr;
	//para.pt_start_infer = nullptr;

    para.data1D = new double[para.dim*para.samplesize]; // new
	para.data1D_vali = new double[para.dim*para.valisize];
	if(p.dim == 1) {
		//para.data1D_infer = new double[para.dim*para.infersize];
	}
    para.pt_start = new int[para.samplesize*para.n];
	para.pt_start_vali = new int[para.valisize*para.n];
	if(p.dim ==1) {
		//para.pt_start_infer = new int[para.infersize*para.n];
	}
    //    bool pathct=false;
    //    bool discrete=true;

    //    double smoothneighbordist=0.001;
    // density to be learned at the current step, linxi 
	//   vector<OnePartition_data>  bestlevelPs;
    OnePartition_data  bestlevelPs;
    //    int bestind;
	
	ensemble.clear();
	
	unsigned seed = time(0);
	srand(seed);
	
	for(int l=0; l< ensemble_sizeP; l++) {
		cout<<"-----------------------------"<<endl;
		cout<<"------------l="<<l<<"--------------"<<endl;
		cout<<"-----------------------------"<<endl;
		// Output related, deleted to make the computation time shorter.
		/*
		string tempdir = ofilename+"RF/tree_"+std::to_string(l)+"/";
		string cline = "mkdir -p "+tempdir;
		const int dir= system(cline.c_str());
		if (dir< 0) return;
		*/
		OneDensity_data onefit;
		onefit.shift.clear();
		onefit.joint_shift.clear();
		vector<int> shift_int(p.dim, 0);
	
		
		//string tempshift = tempdir+"Mshift.txt";
		//ofstream outfileMshift(tempshift.c_str());
		for(int d=0; d <p.dim; d++) {
			shift_int[d] = 0; //rand_int(0, (n_shiftM-1));
			onefit.shift.push_back( (double)shift_int[d] * shift_sizeM );
			//cout<<"For marginal distribution, shift along dimension "<<d<<" is: "<<onefit.shift[d]<<"\t"<<shift_int[d]<<endl;
			//outfileMshift<<onefit.shift[d]<<endl;
		}
		//outfileMshift.close();
	
		//tempshift = tempdir +"joint_shift.txt";
		//ofstream outfilePshift(tempshift.c_str());
		for(int d=0; d<p.dim; d++) {
			shift_int[d] = rand_int(0, (n_shiftP-1));
			onefit.joint_shift.push_back( (double)shift_int[d] * shift_sizeP);
			//cout<<"For joint distribution, shift along dimension "<<d<<" is: "<<onefit.joint_shift[d]<<"\t"<<shift_int[d]<<endl;
			//outfilePshift<<onefit.joint_shift[d]<<endl;
		}
		//outfilePshift.close();
		
		// Output related, deleted to make the computation time shorter.
		//string temptrans = tempdir+"trans_data.txt"; //file for Copula transformed data
		//cout<<temptrans<<endl;
		
		// Output related, deleted to make the computation time shorter.
		onefit.marginalP.clear();
	    //ofstream outfiletrans (temptrans.c_str());
	    for(int d=0; d<dim; d++){
	        cerr<<"dim="<<d<<endl<<endl;
	        para.data.clear();
			para.validata.clear();
	        vector<double> onepiece(1,0); //initialization for a vector of doubles, linxi
	        vector<double> onevariable(para.samplesize, 0);
			vector<double> onevariable_vali(para.valisize, 0);
			//vector<double> oneshift(1,0);
	
	
	        // shift the data along the dth dimension, linxi
			// first shift the original data, 
			// then each coordinate is divided by 2, linxi
			// update both para.data and para.data1D, linxi
			// p.data remains unchanged, linxi
			for(int i=0; i<para.samplesize; i++)   {
	            onepiece[0]=(p.data[i][d]+onefit.shift[d])*0.5;
	            onevariable[i]= (p.data[i][d]+onefit.shift[d])*0.5;
	            para.data.push_back(onepiece);

	            para.data1D[i*para.dim]=(p.data[i][d]+onefit.shift[d])*0.5; // new
	            for(int j=0; j<para.n; j++)
	                para.pt_start[j*para.samplesize+i] = i*para.dim;
	        }
			
	        for(int i=0; i<para.valisize; i++)   {
	            onepiece[0]=(p.validata[i][d]+onefit.shift[d])*0.5;
	            onevariable_vali[i]= (p.validata[i][d]+onefit.shift[d])*0.5;
	            para.validata.push_back(onepiece);

	            para.data1D_vali[i*para.dim]=(p.validata[i][d]+onefit.shift[d])*0.5; // new
	            for(int j=0; j<para.n; j++)
	                para.pt_start_vali[j*para.valisize+i] = i*para.dim;
	        }
	

			// fit BSP for marginal dist, linxi
	        SIS(para, true, bestlevelPs, betaM, Cdepth, 100, "RF");
	
			// print top 5 smallest edge length;
			if(l/50==0) get_smallest_edge_length(bestlevelPs, 5);

	        onefit.marginalP.push_back( bestlevelPs);
			//    print_partition(bestlevelPs);
	        print_partition(onefit.marginalP[d]); // output the marginal dist to terminal, linxi

	        // get marginal quantiles (copula), linxi
	        vector<double> newx = Ftransform(onevariable, onefit.marginalP[d]); // defined in output.h
			// still scale to [0, 0.5] after marginal copula
			scale_by_two(&newx[0], para.samplesize, para.dim, true);
	
			// checkings, linxi
			cerr<<"orginal data is: ";
	        for(int i=0; i<10; i++) cerr<<p.data[i][d]<<" ";
	        cerr<<endl;

	        cerr<<"transformed data is: ";
			for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
	        cerr<<endl;

			// checkings, inverse transformation, linxi
			scale_by_two(&newx[0], para.samplesize, para.dim, false);
	        vector<double> recoveronevariable = inv_Ftransform(newx, onefit.marginalP[d]);
			scale_by_two(&recoveronevariable[0], para.samplesize, para.dim, false);

			cerr<<"recoverd data is: ";
	        for(int i=0; i<10; i++) cerr<<recoveronevariable[i]-onefit.shift[d]<<" ";
	        cerr<<endl;
	
			// update parameter set p for multivariate density estimation
			// updata the data component with the one after copula transformation, linxi
			// output the transformed data
			// transformed data is in [0,1], linxi
	        for(int i=0; i<p.samplesize; i++){
	            //dataF[i][d]=newx[i];            //new
	            //p.data[i][d] = newx[i];       //- by ed520, we use below
	            p.data1D[i*p.dim+d] = newx[i];  //+ by ed520
				// Output related, deleted to make the computation time shorter.
				//if(i< 100) outfiletrans << newx[i]<<'\t';
	        }
			
			vector<double> newx_vali = Ftransform(onevariable_vali, onefit.marginalP[d]);
		
	        for(int i=0; i<p.valisize; i++){
	            //dataF[i][d]=newx[i];            //new
	            //p.data[i][d] = newx[i];       //- by ed520, we use below
	            p.data1D_vali[i*p.dim+d] = newx_vali[i];  //+ by ed520
				// Output related, deleted to make the computation time shorter.
	            //if(i<100) outfiletrans << newx_vali[i]<<'\t';
	        }
	        //outfiletrans<<'\n';
		}
	
		//outfiletrans.close();
	
		if(p.dim > 1) {
			//cout<<"onefit.joint_shift[0]: "<<onefit.joint_shift[0]<<" "<<onefit.joint_shift[1]<<endl;
			//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
			shift_data(p.data1D, p.samplesize, p.dim, onefit.joint_shift);
			//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
			scale_by_two(p.data1D, p.samplesize, p.dim, true);
			
			shift_data(p.data1D_vali, p.valisize, p.dim, onefit.joint_shift);
			scale_by_two(p.data1D_vali, p.valisize, p.dim, true);
			//cout<<p.data1D[0]<<" "<<p.data1D[1]<<" "<<p.data1D[p.dim]<<" "<<p.data1D[p.dim+1]<<endl;
		    //    para.dim = dim;
		    //    para.data.clear();
		    //    para.data=dataF;
		    //    para.levels = p.levels;
		    //    para.n = p.n;
		    //    para.steps = p.steps;
		    //    para.resampling = p.resampling;          //new
		    //   SIS(para, true, bestlevelPs, bestind);

		    // fit the multivariate BSP after copula trans, linxi
			//onefit.jointP.clear();
			SIS(p, true, onefit.jointP, betaP, Cdepth, 100, "RF");

			// print top 5 smallest edge length;
			if(l/100 ==0) get_smallest_edge_length(onefit.jointP, 5);

			// output marginal dist's to files, linxi
		    //     print_Finvtransform(P, marginalP);
			// Output related, deleted to make the computation time shorter.
			/*
		    for (int d = 0; d < dim; d++) {
		        string temp= tempdir+"dim"+std::to_string(d)+".txt";
		        ofstream outfile (temp.c_str());

		        // fprintf(pFile, "d= %d\n", d);
		        print_partition(outfile, onefit.marginalP[d]);
		        outfile.close();
		    }
			*/
			
			// Output related, deleted to make the computation time shorter.
			/*
			string distTrans= tempdir+"transJoint.txt";
			string distOrig=tempdir+"origJoint.txt";	//+ 412
			ofstream outdistTrans(distTrans.c_str());
			ofstream outdistOrig(distOrig.c_str());		//+ 412
			*/
			OnePartition_data &op = onefit.jointP;		

			//    outfile <<endl << "After transform back" << endl;
			for (int i = 0; i < (int) op.sregs.size(); i++) {
				pair<double, double> range;
			    for (int d = 0; d < dim; d++) {
					range = convert_ranges(op.sregs[i].reg_code[d]);
			        double ub,lb;
			        ub = inv_Ftransform(range.second*2-onefit.joint_shift[d], onefit.marginalP[d]);
			        // ub = ub * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
					ub = ub * 2.0 -onefit.shift[d] - 0.000001;
			        lb =  inv_Ftransform(range.first*2-onefit.joint_shift[d], onefit.marginalP[d]);
			        // lb = lb * 2.0 * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
					lb = lb *2.0 - onefit.shift[d] - 0.000001;
			        //         outfile<<range.first<<" "<<range.second<<" ";
					// Output related, deleted to make the computation time shorter.
			        //outdistOrig <<lb << " " << ub <<" ";
			        //outdistTrans << range.first << " " << range.second << " ";
				}
				// calculate the density on a subregion, linxi
				// Output related, deleted to make the computation time shorter.
			    //outdistTrans<< (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;
			  	//outdistOrig << (double)op.sregs[i].num/(double)p.samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)+(double)p.dim*log(2.0)) <<" "<< op.sregs[i].num<<endl;
			}
			// Output related, deleted to make the computation time shorter.
			//outdistTrans.close();
			//outdistOrig.close();	
		}else {
			/*
			string distOrig=tempdir+"origJoint.txt";	//+ 412
			ofstream outdistOrig(distOrig.c_str());	
			print_partition(outdistOrig, onefit.marginalP[0]);
			outdistOrig.close();
			*/	
		}
		// output joint dist to files, linxi	
	
		ensemble.push_back(onefit);			
	}
}
