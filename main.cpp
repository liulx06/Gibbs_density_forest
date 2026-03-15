/*
 Bayesian Sequential Partitioning
 2013-1  by LuoLu
 */

// hello

#include "tree.h"
#include "readdata.h"
#include "examples.h"
#include "timer_ed520.h"    //+ by ed520
#include "density_tree_forest.h"
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
TIMER timer;                //+ by ed520

int testF (vector< map< OnePartition , double, CompairSRegs> >& pathctmap); //don't been called anywhere, thchiu
int BSP_C(vector<string> params);
int BSP_NC(vector<string> params);
int Gibbs_C_shift(vector<string> params);
void print_usage_and_exit();
int density(vector<string> & params);
int density_2(vector<string> & params);     //+ by 412;
int density_tree(vector<string> & params);  //+ by 412;
int BSP_C_forden(vector<string> params);
int BSP_tree_C(vector<string> params);      //+ by 520;
int BSP_tree_NC(vector<string> params);     //+ by 520;
int CDF(vector<string> & params);           //+ by 412;
int CDF_tree(vector<string> & params);	    //+ by 412;
void test_coding(tree BspTree, string ofilename);	//+ by 412; for testing
int dec2bin(int n);					//+ by 412; for testing
void test_partition_generation( Partition &mytestpartition ,  Region_Node *  mytestregion_a    );

int countfrompartition(vector<string> params);

int main(int argc, char** argv) {

 //   cerr << "**" + c::PROG_NAME << '\n';
 //   cerr << "Version: " << c::BUILD << "\n";
    //cerr << "Luo Lu\n";
    //cerr << "http://www.stanford.edu/group/wonglab" << "\n";
    //cerr << "\n";

    string mode = "";
    string temp = "";
    int error_num = 0;
    vector<string> params;

    if (argc < 2) {
        print_usage_and_exit();
    }

    mode = argv[1];
    for (int i = 2; i < argc; i++) {
                temp = argv[i];
                params.push_back(temp);
    }

	//timer.Begin();
	
    if (mode == "BSP_NC") {
         error_num = BSP_NC(params);
    } else if(mode == "BSP_C"){
         error_num= BSP_C_forden(params);
	} else if (mode == "Gibbs_C_shift") {
		 error_num = Gibbs_C_shift(params);
    } else if (mode == "density") {
         error_num = density(params);
    } else if (mode == "density_2") {
         error_num = density_2(params);
    } else if (mode == "density_tree") {
         error_num = density_tree(params);
	} else if(mode == "count"){
        error_num = countfrompartition(params);
    } else if(mode == "BSP_tree_C"){
        error_num=BSP_tree_C(params);
    } else if(mode == "BSP_tree_NC"){
        error_num=BSP_tree_NC(params);
    } else if(mode == "CDF") {
        error_num = CDF(params);
    } else if(mode == "CDF_tree") {
        error_num = CDF_tree(params);
    }
    else{
        print_usage_and_exit();
    }
    //timer.End();
    //timer.RunTimeProfile();

    return error_num;

}


int BSP_C(vector<string> params){
    string usage_text = "Usage: BSP_C <data_file> <output_file> <level1=1000> <level2=200> \n data_file -- One sample each row; MAP partitions output to output_file. Log to STDOUT; level1 -- number of SIS levels for joint distribution; level2 -- number of SIS levels for marginal distribution";

    if (params.size() > 4 || params.size() < 2) {
        cerr << usage_text << endl;
        return 3;
    }
    vector<double> mmax, mmin;
    string ofilename = params[1];
    parameters  p;
    p.data = read_data(params[0], true, mmax, mmin);
    p.dim=(int)p.data[0].size();
 //   p.levels=strTo<int>(params[1]);
    if(params.size()==4) p.levels= strTo<int>(params[2]);
    else p.levels = 1000;
    p.maxpercentage = 0.9;
    p.n = 200;
    p.resampling=2;
    p.samplesize= (int)p.data.size();;
    p.steps=4;
    p.data1D = new double[p.dim*p.samplesize]; //+ by thchiu
    p.pt_start = new int [p.samplesize * p.n]; //+ by thchiu
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }

    int dim = p.dim;
    int samplesize = p.samplesize;
    parameters para;
    para.maxpercentage = 0.9;
    para.samplesize = samplesize;
    para.steps = p.steps;
    para.dim = 1;
    para.n = 200;
    para.data1D = new double[para.dim*para.samplesize]; //+ by thchiu
    para.pt_start = new int [para.samplesize * para.n];//+ by thchiu

    if(params.size()==4)    para.levels = strTo<int>(params[3]);
    else para.levels=200;
    para.resampling = 2;
    //vector<vector<double> >  dataF = p.data;   //commented by thchiu

    vector<OnePartition_data> marginalP;

     //vector<OnePartition_data>  bestlevelPs;  //commented by thchiu
     OnePartition_data bestlevelPs; //added by thchiu
     //int bestind; //commented by thchiu

    for(int d=0; d<dim; d++){
        cout<<"dim="<<d<<endl<<endl;
        para.data.clear();
        vector<double> onepiece(1,0);
        vector<double> onevariable(samplesize, 0);
        for(int i=0; i<samplesize; i++)   {
            onepiece[0]=p.data[i][d];
            onevariable[i]= p.data[i][d];
            para.data.push_back(onepiece);
            para.data1D[i*para.dim]=p.data[i][d];   //+ by thchiu
            for(int j=0; j<para.n; j++)             //+ by thchiu
                para.pt_start[j*para.samplesize+i] = i*para.dim;
        }


      
        SIS(para, true, bestlevelPs);   //thchiu
        marginalP.push_back( bestlevelPs);//thchiu
        print_partition(marginalP[d]);

  //      print_partition(marginalP[d]);
        vector<double> newx = Ftransform(onevariable, marginalP[d]);

        for(int i=0; i<10; i++) cerr<<onevariable[i]<<" ";
        cerr<<endl;

        for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
        cerr<<endl;

        vector<double> recoveronevariable = inv_Ftransform(newx, marginalP[d]);

        for(int i=0; i<10; i++) cerr<<recoveronevariable[i]<<" ";
        cerr<<endl;

        //for(int i=0; i<samplesize; i++)        //dataF[i][d]=newx[i];   //commented by thchiu
        //    p.data[i][d]=newx[i];

        for(int i=0; i<samplesize; i++){        //+ by thchiu
            p.data1D[i*p.dim+d] = newx[i];
        }
    }
    

     //SIS(para, true, bestlevelPs, bestind);   //commented by thchiu
     SIS(p, true, bestlevelPs);    //added by thchiu

     ///////// Output file is commented by ed520  ////////////
#define OUTPUT
#ifdef OUTPUT
    //     print_Finvtransform(P, marginalP);
    for (int d = 0; d < dim; d++) {
        string temp= ofilename+'_'+toStr<int>(d)+".txt";
        ofstream outfile (temp.c_str());
        // fprintf(pFile, "d= %d\n", d);
        print_partition(outfile, marginalP[d]);
        outfile.close();
    }
    string temp= ofilename+"_Big.txt";
    ofstream outfile (temp.c_str());

    //OnePartition_data op = bestlevelPs[bestind];  //commented by thchiu
    OnePartition_data &op = bestlevelPs; //added by thchiu

     //print_partition(outfile, bestlevelPs[bestind]);  //commented by thchiu
    print_partition(outfile, bestlevelPs);  //added by thchiu

    outfile <<endl << "After transform back" << endl;
    for (int i = 0; i < (int) (int) op.sregs.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < dim; d++) {
            range = convert_ranges(op.sregs[i].reg_code[d]);
            double ub,lb;
            ub = inv_Ftransform(range.second, marginalP[d]);
            ub = ub * (mmax[d] - mmin[d])+ mmin[d];
            lb =  inv_Ftransform(range.first, marginalP[d]);
            lb = lb * (mmax[d] - mmin[d])+ mmin[d];
            outfile <<lb << " " << ub <<" ";
        }
        outfile<< op.sregs[i].num<<endl;

    }

    outfile.close();

 // output the overall regions without copula

    temp= ofilename+"_Big_total.txt";
    ofstream outfile2 (temp.c_str());
    for (int i = 0; i < (int) (int) op.sregs.size(); i++) {
          pair<double, double> range;
          vector<usint_mask> reg_codes= op.sregs[i].reg_code;
          vector<pair<double, double> > regs;
          for (int d = 0; d < dim; d++) {

              range = convert_ranges(reg_codes[d]);
              pair<double, double> bounds;
              bounds.first = inv_Ftransform(range.first, marginalP[d]);
              bounds.second =  inv_Ftransform(range.second, marginalP[d]);
              regs.push_back(bounds);
          }

          recover_fullregions(regs, marginalP,mmax, mmin, outfile2);

    }

    outfile2.close();
#endif

    return 0;

}


// output the partition after normalization
int BSP_NC(vector<string> params){

//    string usage_text = "Usage: BSP_NC <data_file> <output_file> <level1=1000>\n"
 //           + " data_file -- One sample each row\n"
//            + " MAP partitions output to output_file. Log to STDOUT \n"
//            + " level -- number of SIS levels for joint distribution";


    if (params.size() > 4 || params.size()<2) {
//        cerr << usage_text << endl;
        cerr<<"wrong input!\n";
        return 3;
    }
    string ofilename = params[1];
    parameters p;
    vector<double> mmax, mmin;
    p.data = read_data(params[0], true, mmax, mmin);
    cerr<<p.data[0][0]<<" "<<p.data[1][0]<<" "<<p.data[2][0]<<" "<<p.data[3][0]<<endl;

    p.dim=(int)p.data[0].size();
    if(params.size()>2)   p.levels=strTo<int>(params[2]);
    else p.levels = 1000;
    p.maxpercentage = 0.9;
    if(params.size()>3)   p.n=strTo<int>(params[3]);   //+ by ed520
    else p.n = 200;
    p.resampling=2;
    p.samplesize= (int)p.data.size();;
    p.steps=4;

    p.data1D = new double[p.dim*p.samplesize];  // new

    for(int d = 0;d<p.dim;d++){
        for(int i = 0;i<p.samplesize;i++){
            p.data1D[i*p.dim+d] = p.data[i][d];

        }
    }


    p.pt_start = new int [p.samplesize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }

  
    OnePartition_data Ps;   //added by thchiu

    //SIS(p, true, Ps, bestind);    //commented by thchiu
    SIS(p, true, Ps, ofilename, mmax, mmin);   //added by thchiu
    
//    print_partition(Ps);
    //ofstream outfile (ofilename.c_str()); //- by ed520
    string temp= ofilename+"_Big.txt";  //+ by ed520
    ofstream outfile (temp.c_str()); //+ by ed520
 //   print_partition(outfile,Ps[bestind]);

  //  outfile <<endl<<endl;


    //print_partition(outfile,Ps[bestind], mmax, mmin); //commented by thchiu
    print_partition(outfile,Ps, mmax, mmin);    //added by thchiu

    outfile.close();
    //+ by ed520
    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<p.dim<<endl;
    cout<<"Data_num = "<<p.samplesize<<endl;
    cout<<"# of Partition = "<<p.n<<endl;
    //+ end by ed520
    return 0;

}


int density(vector<string> & params){

 //  string usage_text = "Usage: \n"
 //           + " density <-c/-n> <Partitions> <sample_data> \n"
 //           + "       -c/-n -- c=use copula, n=no copula\n"
 //           + "  partitions -- partitions of a distribution\n"
 //           + " sample_data -- Each row one data point\n "
 //           + " Output the density.\n";

    if (params.size() != 4) {
        cerr<<params.size()<<endl;
        cerr << "usage_text" << endl;
        return 3;
    }

//    cerr<<params[1]<<'\n'<<params[2]<<''
    vector<vector<double> > test_data = read_data(params[1]);  //data
    vector<vector<double> > mmax_min = read_data(params[3]); // file size: dim*2

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();
    vector<double> dens(test_N, 0);

    for(int i=0; i<test_N; i++){
        for(int d=0; d<dim; d++){
       test_data[i][d] = (test_data[i][d] - mmax_min[d][1]) / (mmax_min[d][0] - mmax_min[d][1]);
        }
    }
    cerr << test_N << " data points in " << dim << " dimensions.\n";
    cerr << "some data: "<<test_data[10][1]<<" "<<test_data[100][20]<<endl; //
    if (params[0]=="-c"){
         vector<vector<vector<double> > > MarginalPs;
         for(int d=0; d<dim; d++){
               MarginalPs.push_back(read_data( params[2]+'_'+toStr<int>(d)+".txt"));
         }
         vector<vector<double> > P = read_data( params[2]+"_Big.txt");//after copula transback, after normalize
         cerr<<"P="<<P[0][2]<<" "<<P[1][5]<<endl;

         for (int i = 0; i < test_N; i++) {
             if(max(test_data[i])>=1 || min(test_data[i])<=0){ // outside cubic
               //  cout<<i<<'\t'<<dens[i]<<'\n';
                 cerr<<"0"<<endl; //
                 continue;
             }
             dens[i] = densities_from_partition(test_data[i], P);  //the 2*dim col of P should be density

            if(i<10) cerr<<"\n i= \t"<<dens[i]<<endl;
            for (int d = 0; d < dim; d++) {
                vector<double> onepiece;
                onepiece.push_back(test_data[i][d]);
                dens[i] *= densities_from_partition(onepiece, MarginalPs[d]);
                 if(i<10) cerr<<dens[i]<<'\t';
            }

        }
    }
    else  if (params[0]=="-n"){
        vector<vector<double> > P = read_data( params[1]);
         for (int i = 0; i < test_N; i++) {
            dens[i] = densities_from_partition(test_data[i], P);
        }
    }
    else {
        cerr << "usage_text" << endl;
        return 1;
    }

    string tempout= params[2]+"_den.txt";
    ofstream outfileout (tempout.c_str());
    for(int i=0; i<test_N; i++){
        for(int d=0; d<dim; d++){
             dens[i]/= (mmax_min[d][0] - mmax_min[d][1]);

        }
        outfileout <<dens[i]<<'\n';
    }

    outfileout.close();


    return 0;

}


int density_2(vector<string> & params){

    if (params.size() != 4) {
        cerr<<params.size()<<endl;
        cerr << "usage_text" << endl;
        return 3;
    }

    vector<vector<double> > test_data = read_data(params[1]);  //data

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();
    vector<double> dens(test_N, 0);

    cerr << test_N << " data points in " << dim << " dimensions.\n";
    
    if (params[0]=="-c"){
         vector<vector<vector<double> > > MarginalPs;
         vector<vector<double> > mmax_min = read_data( params[2] + "_maxmin.txt");
         
         for (int i=0; i<test_N; i++){
             for (int d=0; d<dim; d++){
                 test_data[i][d] = (test_data[i][d] - mmax_min[d][1]) / (mmax_min[d][0] - mmax_min[d][1]);
             }
         }
         for(int d=0; d<dim; d++){
               MarginalPs.push_back(read_data( params[2]+'_'+toStr<int>(d)+".txt"));
         }
                                                                          
         vector<vector<double> > P = read_data( params[2]+"_Big.txt");//after copula transback, after normalize

         for (int i = 0; i < test_N; i++) {
            if(max(test_data[i])>=1 || min(test_data[i])<=0){ // outside cubic
                cerr<<"0: This data point is outside the cube."<<endl; //
                continue;
            }
            dens[i] = densities_from_partition(test_data[i], P);  //the 2*dim col of P should be density

            for (int d = 0; d < dim; d++) {
                vector<double> onepiece;
                onepiece.push_back(test_data[i][d]);
                dens[i] *= densities_from_partition(onepiece, MarginalPs[d]);
            }
            for (int d=0; d<dim; d++){
                dens[i] /= (mmax_min[d][0] - mmax_min[d][1]);                
            }
        }
    }
    else  if (params[0]=="-n"){
        vector<vector<double> > P = read_data( params[2]+"_Big.txt");
         for (int i = 0; i < test_N; i++) {
            dens[i] = densities_from_partition(test_data[i], P);
        }
    }
    else {
        cerr << "usage_text" << endl;
        return 1;
    }

    string tempout= params[3]+"_den.txt";       // modified 2 -> 3
    ofstream outfileout (tempout.c_str());
    for(int i=0; i<test_N; i++){
        outfileout <<dens[i]<<'\n';
    }

    outfileout.close();

    return 0;

}


int density_tree(vector<string> & params){

    if (params.size() != 4) {
        cerr<<params.size()<<endl;
        cerr << "usage_text" << endl;
        return 3;
    }

    vector<vector<double> > test_data = read_data(params[1]);  //data

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();
    vector<double> dens(test_N, 0);

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    if (params[0]=="-c"){
        vector<vector<vector<double> > > MarginalPs;
        vector<vector<vector<double> > > Tree_marginals;
        
        vector<vector<double> > mmax_min = read_data( params[2] + "_maxmin.txt");
        vector<vector<double> > Fx = test_data;	//+ 412
        for (int i=0; i<test_N; i++){
            for (int d=0; d<dim; d++){
                test_data[i][d] = (test_data[i][d] - mmax_min[d][1]) / (mmax_min[d][0] - mmax_min[d][1]);
            }
        }
        for(int d=0; d<dim; d++){
            MarginalPs.push_back(read_data( params[2]+'_'+toStr<int>(d)+".txt"));
            Tree_marginals.push_back( read_data( params[2]+'_'+toStr<int>(d)+"_Tree.txt"));
        }
        vector<vector<double> > P_copula = read_data( params[2]+"_Tree.txt");
                                                                                                                               
        for (int i = 0; i < test_N; i++) {
            dens[i] = 1;
            for (int d=0; d<dim; d++){
                vector<double> onepiece;
                onepiece.push_back(test_data[i][d]);
                Fx[i][d] = CDF_from_partition(onepiece, MarginalPs[d]);
                
                dens[i] *= densities_from_partition(onepiece, MarginalPs[d]);	// 1D Marginal PDF, look up region table is quite simple
                //dens[i] *= densities_from_tree(onepiece, Tree_marginals[d]);	// Using Marginal PDFs' tree files is OK, but so far these files are not ready
            }

            if(max(test_data[i])>=1 || min(test_data[i])<=0){ // outside cubic
                cerr<<"0: This data point is outside the cube. "<<endl; //
                continue;
            }
            
            dens[i] *= densities_from_tree(Fx[i], P_copula);
            //cout << densities_from_tree(Fx[i], P_copula) <<'\t';
            
            for (int d=0; d<dim; d++){
                dens[i] /= (mmax_min[d][0] - mmax_min[d][1]);                
            }
        }
    }
    else  if (params[0]=="-n"){
        vector<vector<double> > P = read_data( params[2]+"_Tree.txt");
        
        for (int i = 0; i < test_N; i++) {
            dens[i] = densities_from_tree(test_data[i], P);
        }
    }
    else {
        cerr << "usage_text" << endl;
        return 1;
    }
    string tempout= params[3]+"_den.txt";       
    ofstream outfileout (tempout.c_str());
    for(int i=0; i<test_N; i++){
        outfileout <<dens[i]<<'\n';
    }

    outfileout.close();
    
    return 0;

}


int countfrompartition(vector<string> params){

//    string usage_text = "Usage: count dim <data_file> <partition_file><output number file> \n"
 //           + " data_file -- One sample each row\n";


    if (params.size() > 4 || params.size()<4) {
//        cerr << usage_text << endl;
        cerr<<"wrong input!\n";
        return 3;
    }
    int dim = atoi(params[0].c_str());
    cerr<<"dim="<<dim<<endl;
    string datafilename = params[1];
    string partitionfilename = params[2];
    string ofilename = params[3];
    string line;
    ifstream infile;

    infile.open(partitionfilename.c_str());
    vector<vector<double> > partition;

    while (!infile.eof()) {
        getline(infile, line);
        trim2(line);
        if (line.length() == 0) continue;
        vector<string> ll = split(line);
        vector<double> d;
        d.resize(2*dim);
        for (int i = 0; i < 2*dim; i++) {
            d[i] = strTo<double>(ll[i]);
        }
  //      cerr<<d[7]<<" ";
        partition.push_back(d);
    }
//    cerr<<endl<<partition[3][7]<<endl;
    infile.close();
    int regnum = (int)partition.size();
    infile.open(datafilename.c_str());
    vector<int> nums(regnum, 0);

    while (!infile.eof()) {
        getline(infile, line);
        trim2(line);

        if (line.length() == 0) continue;

        vector<string> ll = split(line);

        if ((int) ll.size() != dim) {
            cerr << "ERROR: Bad line dim: " << line << '\n';
            exit(2);
        }

        vector<double> d;
        d.resize(dim);

        for (int i = 0; i < dim; i++) {
            d[i] = strTo<double>(ll[i]);
        }
        for(int m=0; m<regnum; m++){
            bool inreg=true;
            for(int j=0; j<dim; j++){
                if(d[j]> partition[m][2*j+1] || d[j]<= partition[m][2*j]){
                    inreg=false;
                    break;
                }
            }
            if(inreg) {
                nums[m]++;
                break;
            }

        }
    }
    infile.close();

    ofstream outfile (ofilename.c_str());

    for(int i=0; i<regnum; i++) outfile<<nums[i]<<" ";
    outfile<<endl;
    cerr<< endl << sum(nums) <<endl;
    return 0;
}


int BSP_C_forden(vector<string> params){
    string usage_text = "Usage: BSP_C <data_file> <output_file> <level=1000> <number of partitions=200> \n data_file -- One sample each row; MAP partitions output to output_file. Log to STDOUT; level -- number of SIS levels for marginal/joint distributions; number of partitions -- number of partitions to sample from the posterior";

	// changed the length of input, linxi
    if (params.size() > 5 || params.size() < 2) {
        cerr << usage_text << endl;
        return 3;
    }
    string ofilename = params[1];
    // the parameter set for joint dist, linxi
	parameters  p;
    vector<double> mmax, mmin;
	
	// read in data, stored in p.data, linxi
	// Data will be scaled along each direction, linxi
	// ranges of data are stored in mmax, mnin, the two parameters will be changed after the function is called, linxi 
    // p.data = read_data(params[0], true, mmax, mmin);
	
	// directly generate data calling functions, linxi
	int datasize;
	datasize = strTo<int>(params[0]);
	sample_normal(datasize, p.data);
	scale_data(p.data, 0.0, 1.0, true, false, mmax, mmin);
	//mmax.resize(2);
	//mmin.resize(2);
	//mmax[0] = 1.0;
	//mmax[1] = 1.0;
	//mmin[0] = 0.0;
	//mmin[1] = 0.0; 
	
	// this is test data for calculating KL dist, linxi
	int testsize;
	testsize = strTo<int>(params[4]);
	sample_normal(testsize, p.testdata);
 
	// output the range of data, linxi
    string tempmaxmin= ofilename+"_maxmin.txt";
    ofstream outfilemaxmin (tempmaxmin.c_str());
    for(int d=0; d<(int)p.data[0].size();d++){
        outfilemaxmin << mmax[d]<<'\t'<<mmin[d]<<'\n';
    }
    outfilemaxmin.close();

    p.dim=(int)p.data[0].size();
 	// p.levels=strTo<int>(params[1]);
	// MAXIMUM number of subregions in a partition cannot be larger than 1000, linxi
    if(params.size()==4) p.levels= strTo<int>(params[2]);
    else p.levels = 1000;
	// parameter for smoothing, the maximum probability for a potential cut cannot exceed 0.9, linxi
	// in order to maintain diversity of partitions, linxi
    p.maxpercentage = 0.9;
	// MAXIMUM number of partitions in an ensemble cannot be larger than 200, linxi
    if(params.size()==4) p.n = strTo<int>(params[3]); //+ by ed520
    else p.n = 200;
	// flag for turning on resample(), linxi
    p.resampling=2;
	// sample size, linxi
    p.samplesize= (int)p.data.size();;
	// resampling every "assigned (eg. 4 here)" steps, linxi
    p.steps=4;

    p.data1D = new double[p.dim*p.samplesize];  // new
    p.pt_start = new int [p.samplesize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }

    //    bool pathct=false;
    //    bool discrete=true;

    //    double smoothneighbordist=0.001;

    int dim = p.dim;
    int samplesize = p.samplesize;
    parameters para; // the papameter set for marginal dist, linxi
    para.maxpercentage = 0.9;
    para.samplesize = samplesize;
    para.steps = p.steps;
    para.dim = 1;
    //para.n = 200; //- by ed520
    para.n = p.n; //+ by ed520
    para.data1D = new double[para.dim*para.samplesize]; // new
    para.pt_start = new int[para.samplesize*para.n];

    //if(params.size()==4)    para.levels = strTo<int>(params[3]);  //- by ed520
    //else para.levels=200; //- by ed520
    para.levels = p.levels; //+ by ed520
    para.resampling = 2;
    //    vector<vector<double> >  dataF= p.data;
	// partitions for marginal dist's, linxi
    vector<OnePartition_data> marginalP; // struct defined in sampling.h, linxi

    // density to be learned at the current step, linxi 
	//   vector<OnePartition_data>  bestlevelPs;
    OnePartition_data  bestlevelPs;
    //    int bestind;
    
	string temptrans = ofilename+"_trans.txt"; //file for Copula transformed data
    ofstream outfiletrans (temptrans.c_str());
    for(int d=0; d<dim; d++){
        cerr<<"dim="<<d<<endl<<endl;
        para.data.clear();
        vector<double> onepiece(1,0); //initialization for a vector of doubles, linxi
        vector<double> onevariable(samplesize, 0);
		
		
        for(int i=0; i<samplesize; i++)   {
            onepiece[0]=p.data[i][d];
            onevariable[i]= p.data[i][d];
            para.data.push_back(onepiece);

            para.data1D[i*para.dim]=p.data[i][d]; // new
            for(int j=0; j<para.n; j++)
                para.pt_start[j*para.samplesize+i] = i*para.dim;
        }
		

		// fit BSP for marginal dist, linxi
        SIS(para, true, bestlevelPs);
		
		// print top 5 smallest edge length;
		get_smallest_edge_length(bestlevelPs, 5);

        marginalP.push_back( bestlevelPs);
		//    print_partition(bestlevelPs);
        print_partition(marginalP[d]); // output the marginal dist to terminal, linxi

        // get marginal quantiles (copula), linxi
        vector<double> newx = Ftransform(onevariable, marginalP[d]); // defined in output.h
		
		// checkings, linxi
		cerr<<"orginal data is: ";
        for(int i=0; i<10; i++) cerr<<onevariable[i]<<" ";
        cerr<<endl;

        cerr<<"transformed data is: ";
		for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
        cerr<<endl;

		// checkings, inverse transformation, linxi
        vector<double> recoveronevariable = inv_Ftransform(newx, marginalP[d]);

		cerr<<"recoverd data is: ";
        for(int i=0; i<10; i++) cerr<<recoveronevariable[i]<<" ";
        cerr<<endl;
		
		// update parameter set p for multivariate density estimation
		// updata the data component with the one after copula transformation, linxi
        for(int i=0; i<samplesize; i++){
            //dataF[i][d]=newx[i];            //new
            //p.data[i][d] = newx[i];       //- by ed520, we use below
            p.data1D[i*p.dim+d] = newx[i];  //+ by ed520
            outfiletrans << newx[i]<<'\t';
        }
        outfiletrans<<'\n';
    }

    outfiletrans.close();
    //   SIS(para, true, bestlevelPs, bestind);

    // fit the multivariate BSP after copula trans, linxi
	SIS(p, true, bestlevelPs);
	
	// print top 5 smallest edge length;
	get_smallest_edge_length(bestlevelPs, 5);

	// output marginal dist's to files, linxi
    //     print_Finvtransform(P, marginalP);
    for (int d = 0; d < dim; d++) {
        string temp= ofilename+'_'+toStr<int>(d)+".txt";
        ofstream outfile (temp.c_str());

        // fprintf(pFile, "d= %d\n", d);
        print_partition(outfile, marginalP[d]);
        outfile.close();
    }
	
	// output joint dist to files, linxi
    string temp= ofilename+"_Big.txt";
    string temp2=ofilename+"_OriBig.txt";	//+ 412
    ofstream outfile (temp.c_str());
    ofstream outfile2(temp2.c_str());		//+ 412
    //   print_partition(outfile, bestlevelPs[bestind]);
    //   OnePartition_data op = bestlevelPs[bestind];
    OnePartition_data &op = bestlevelPs;

    //    outfile <<endl << "After transform back" << endl;
    for (int i = 0; i < (int) (int) op.sregs.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < dim; d++) {
            range = convert_ranges(op.sregs[i].reg_code[d]);
            double ub,lb;
            ub = inv_Ftransform(range.second, marginalP[d]);
            //     ub = ub * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
            lb =  inv_Ftransform(range.first, marginalP[d]);
            //     lb = lb * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
            //         outfile<<range.first<<" "<<range.second<<" ";
            outfile <<lb << " " << ub <<" ";
            outfile2 << range.first << " " << range.second << " ";
        }
		// calculate the density on a subregion, linxi
        outfile<< (double)op.sregs[i].num/(double)samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;
        outfile2 << (double)op.sregs[i].num/(double)samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;
    }

    outfile.close();
    outfile2.close();
    //+ by ed520
    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<p.dim<<endl;
    cout<<"Data_num = "<<p.samplesize<<endl;
    cout<<"# of Partition = "<<p.n<<endl;
    //+ end by ed520
	
	// calculating KL divergence, linxi
	double KLdist;
	KLdist = KL_fromf(sample_normal_density, p.testdata, op, marginalP, 0.0, 0.0);
	cout<<"KL divergence from the true is: "<<KLdist<<endl;
    return 0;

}


int Gibbs_C_shift(vector<string> params){    
	string usage_text = "Usage: Gibbs_C_shift Please check the ReadMe file";

	// changed the length of input, linxi
	// Include beta as part of input, linxi, 9/8/24
	// Include shape constraint depth as part of input, linxi, 9/9/24
	int total_para = 19; 
    if (params.size() > total_para || params.size() < 2) {
        cerr << usage_text << endl;
        return 3;
    }
    string ofilename = params[1];
    // the parameter set for joint dist, linxi
	parameters  p;
	p.data.clear();
	p.validata.clear();
	p.testdata.clear();
	
	int datasize;
	datasize = strTo<int>(params[0]);
	
	int testsize;
	testsize = strTo<int>(params[4]);
	
	// this is the upper bound for the shift
	// should be of the form 2^(- grid_scale)
	// input is an integer, only the exponent, linxi
	int grid_scaleM;
	if(params.size() == total_para) grid_scaleM = strTo<int>(params[5]);
	else grid_scaleM = 5;
	double grid_sizeM = exp(-(double)(grid_scaleM)*log(2.0));
	double tempNum = exp((double)(grid_scaleM)*log(2.0));
	int n_shiftM = (int)tempNum;
	//cout<<tempNum<<"\t"<<(double)n_shift +1 - tempNum <<"\t"<<tempNum - (double)n_shift<<endl;
	if((double)n_shiftM +1 - tempNum < tempNum - (double)n_shiftM) n_shiftM +=1;
	cout<<"The largest possible shift along each dimension is: "<<grid_sizeM<<"; The number of shifts for learning the marginal distributions is: "<<n_shiftM<<endl;
	double shift_sizeM = grid_sizeM * grid_sizeM;
	cout<<"The unit of shift for learning marginal distributions is: "<<shift_sizeM<<endl;
	
	// shifting unit for joint distribution after the copula transformation
	// This is usually set to be zero in the implementation.
	int grid_scaleP = 0;
	double grid_sizeP = exp(-(double)(grid_scaleP)*log(2.0));
	tempNum = exp((double)(grid_scaleP)*log(2.0));
	int n_shiftP = (int)tempNum;
	//cout<<tempNum<<"\t"<<(double)n_shift +1 - tempNum <<"\t"<<tempNum - (double)n_shift<<endl;
	if((double)n_shiftP +1 - tempNum < tempNum - (double)n_shiftP) n_shiftP +=1;
	cout<<"The largest possible shift along each dimension is: "<<grid_sizeP<<"; The number of shifts for learning the joint distribution is: "<<n_shiftP<<endl;
	double shift_sizeP = grid_sizeP * grid_sizeP;
	cout<<"The unit of shift for learning marginal distributions is: "<<shift_sizeP<<endl;
	
	// Define the sample size for inference 
	double inferP;
	if(params.size()== total_para)
		inferP =strTo<double>(params[6]);
	else
		inferP =1.0; 
	
	// proportion of data for leaf node estimation given the forest topology
	double splitP;
	if(params.size()== total_para)
		splitP =strTo<double>(params[7]);
	else
		splitP =0.5;
	
	//int ensemble_sizeM;
	//if(params.size() ==15) ensemble_sizeM = strTo<int>(params[7]);
	//else ensemble_sizeM = 100;
	
	// The number of samples of GDFs
	int iter_start;
	if(params.size()== total_para) iter_start = strTo<int>(params[8]);
	else iter_start =0;
	
	// total number of shifts, also the number of trees within each forests
	// for high-dimensional cases, if the number of possible shifts is too large,
	// then random pick shifts from all possibles ones, linxi
	int ensemble_sizeP;
	if(params.size() == total_para) ensemble_sizeP = strTo<int>(params[9]);
	else ensemble_sizeP = 20;
	
	// p is parameter for multivariate density	
    // p.dim=(int)p.data[0].size();
	// MAXIMUM number of subregions in a partition cannot be larger than 1000, linxi
    if(params.size()== total_para) p.levels= strTo<int>(params[2]);
    else p.levels = 1000;
	cout<<"number of levels is: "<<p.levels<<endl;
	// parameter for smoothing, the maximum probability for a potential cut cannot exceed 0.9, linxi
	// in order to maintain the diversity of the tree structure, linxi
    p.maxpercentage = 0.9;
	// MAXIMUM number of tree samples cannot be larger than 200, linxi
	// Note: we use particle filters to grow trees under the OPT prior
	// We draw multiple trees samples from the posterior, then draw one sample according to the weights
    if(params.size()== total_para) p.n = strTo<int>(params[3]); //+ by ed520
    else p.n = 200;
	// flag for turning on resample(), linxi
    p.resampling=2;
	// sample size for leaf node estimation, linxi
	if((int)p.data.size()>0) {
		p.samplesize = (int)p.data.size();
	} else {
		p.samplesize= (int)(datasize*splitP);
	} 
	// sample size for risk evaluation, linxi 
	if((int)p.validata.size()>0) {
		p.valisize = (int)p.validata.size();
	} else {
		p.valisize = datasize - p.samplesize;
	}
	
	// Note: sample size for inference should be smaller
	double infersize = (double)(datasize);
	//double scale = log(infersize);
	//infersize = infersize /scale;
	infersize = infersize * inferP; 	
	p.infersize = (int)infersize;
	//p.infersize = p.samplesize + p.valisize;
	cout<<"Sample size is: "<<datasize<<"; training size is :"<<p.samplesize<<"; validation size is: "<<p.valisize<<"; inference dataset is of size: "<<p.infersize<<endl;
	// resampling every "assigned (eg. 4 here)" steps, linxi
    p.steps=4;
	
	// Initialize the pointers, although they will be reallocated later, linxi
	p.dim = 1;
    p.data1D= new double[p.dim*p.samplesize];  // new
	p.data1D_vali = new double[p.dim*p.valisize];
    p.pt_start = new int[p.samplesize * p.n];
	p.pt_start_vali = new int[p.valisize*p.n];
	/*
	// This part can be skipped, will be initialized later, linxi
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
        for(int j=0; j<p.valisize; j++){
            p.pt_start_vali[i*p.valisize+j] = j*p.dim;
        }
    }
	*/
		
    //int dim = p.dim;
    // int samplesize = p.samplesize;
    parameters para; // the papameter set for marginal dist, linxi
    para.maxpercentage = 0.9;
    para.samplesize = p.samplesize;
	para.valisize = p.valisize;
	para.infersize = p.infersize;
    para.steps = p.steps;
    para.dim = 1;
    //para.n = 200; //- by ed520
    para.n = p.n; //+ by ed520
	
	// Initialize the pointers, although they will be reallocated later, linxi
    para.data1D = new double[para.dim*para.samplesize]; // new
    para.pt_start= new int[para.samplesize*para.n];
	para.data1D_vali = new double[para.dim*para.valisize]; // new
    para.pt_start_vali= new int[para.valisize*para.n];
	if(p.dim == 1) {
		para.data1D_infer = new double[para.infersize*para.n];
		para.pt_start_infer = new int[para.infersize *para.n];		   	
	}
	
    //if(params.size()==4)    para.levels = strTo<int>(params[3]);  //- by ed520
    //else para.levels=200; //- by ed520
    para.levels = p.levels; //+ by ed520
    para.resampling = 2;
    //    vector<vector<double> >  dataF= p.data;
	
	int nrep= strTo<int>(params[13]);
	cout<<"Number of repetitions of the experiment is: "<<nrep<<endl;
	
	int RFdepth= strTo<int>(params[14]);
	cout<<"The depth limit for random forests (density) is: "<<RFdepth<<endl;
	
	double betaM;
	if(params.size()==total_para)
		betaM =strTo<double>(params[11]);
	else
		betaM = 1.0;
	
	double betaP;
	if(params.size()==total_para)
		betaP =strTo<double>(params[12]);
	else
		betaP = 1.0;
	
	// The method for density learning
	string learning_method = params[15];
	
	// Whether or not to perform UQ
	bool UQ_flag = (params[16] == "true" || params[16] == "1");
	
	// temperature for sampling trees
	double tree_temp =strTo<double>(params[17]);
	// temperature for sampling forests
	double forest_temp= strTo<double>(params[18]);
	
	int Cdepth = 30; //coupling depth, not used in the current version of code
	
	// Simulate test data for calculation of KL, only generate once, linxi
	string data_model = params[10];
	unsigned seed = 283+testsize/100;
	srand(seed);
	sample_data(testsize, p.testdata, data_model);
	p.dim = (int)p.testdata[0].size();

	string cline = "mkdir -p "+ofilename;
	const int dir= system(cline.c_str());
	if (dir< 0) return 5;
	
	// only keep the test data once
	if(learning_method == "bma" && grid_scaleM == 5) {
		string testdatafile = ofilename+"testdata.txt";
		ofstream outtestdata(testdatafile.c_str());
		for(int i=0; i<testsize; i++) {
			for(int d=0; d<(int)p.testdata[0].size(); d++) {
				outtestdata<<p.testdata[i][d]<<"\t";
			}
			outtestdata<<endl;
		}
		outtestdata.close();	
	}
	
	string KL_BSP = ofilename+"KL_BSP.txt";
	ofstream outKL_BSP(KL_BSP.c_str());

	string KL_BSP_shift = ofilename+"KL_BSP_shift.txt";
	ofstream outKL_BSP_shift(KL_BSP_shift.c_str());

	string KL_BSP_fof = ofilename+"KL_BSP_fof.txt";
	ofstream outKL_BSP_fof(KL_BSP_fof.c_str());

	string KL_BSP_Gibbs = ofilename+"KL_BSP_Gibbs.txt";
	ofstream outKL_BSP_Gibbs(KL_BSP_Gibbs.c_str());
	
	string time_BSP = ofilename+"running_time_BSP.txt";
	ofstream outTime_BSP(time_BSP.c_str());

	string time_BSP_shift = ofilename+"running_time_BSP_shift.txt";
	ofstream outTime_BSP_shift(time_BSP_shift.c_str());

	string time_BSP_fof = ofilename+"running_time_BSP_fof.txt";
	ofstream outTime_BSP_fof(time_BSP_fof.c_str());
	
	string time_UQ= ofilename+"running_time_UQ.txt";
	ofstream outTime_UQ(time_UQ.c_str());
	
	string KL_RF = ofilename+"KL_RF.txt";
	ofstream outKL_RF(KL_RF.c_str());

	string time_RF = ofilename+"running_time_RF.txt";
	ofstream outTime_RF(time_RF.c_str());			
	
	// repeat multiple times to calculate the average KL for each method
	for(int iter=iter_start; iter <iter_start+nrep; iter++) {
		p.data.clear();
		p.validata.clear();
		p.inferdata.clear();
		para.data.clear();
		para.validata.clear();
		para.testdata.clear();
		para.inferdata.clear();
		// read in data, stored in p.data, linxi
		// Data will be scaled along each direction, linxi
		// ranges of data are stored in mmax, mnin, the two parameters will be changed after the function is called, linxi 
	    // // p.data = read_data(params[0], true, mmax, mmin); this is old version before separating scale from reading in files, linxi
		// p.data = read_data_only (params[0]);
		// // Data will be shifted to [0, 0.5]^p for "shift" option
		// scale_data(p.data, 0.0, 1, true, false, mmax, mmin);
		// scale_by_two(p.data, true);	
		
		// Call a function to generate data, added by linxi, Nov 2, 2025
		// "data_model" is the distribution from the input
		seed = iter*17+23+datasize/100;
		srand(seed);
		sample_data(p.samplesize, p.data, data_model);
		sample_data(p.valisize, p.validata, data_model);
		if(p.dim == 1 && UQ_flag) {
			sample_data(p.infersize*ensemble_sizeP, p.inferdata, data_model); // size of data for UQ is the same as that for learning the posterior 
		}
	
		// for(int i=0; i<10;i++) {
		// 	cout<<"One sample is: "<<p.data[i][0]<<" "<<p.data[i][1]<<endl;
		//	cout<<"The density is: "<<sample_normal_AR1_density(p.data[i], 0.6)<<endl;
		// }
		
		string ofilenamerep = ofilename + "rep_"+std::to_string(iter)+"/";
		string cline = "mkdir -p "+ofilenamerep;
		const int dir= system(cline.c_str());
		if (dir< 0) return 5;
		
		// A vector for the range of data
		// No need to update this for simulated data
		vector<double> mmax((int)p.data[0].size(), 1.0), mmin((int)p.data[0].size(), 0.0);

		// output the range of data, linxi
	    string tempmaxmin= ofilenamerep+"range.txt";
	    ofstream outfilemaxmin (tempmaxmin.c_str());
	    for(int d=0; d<(int)p.data[0].size();d++){
	        outfilemaxmin << mmax[d]<<'\t'<<mmin[d]<<'\n';
	    }
	    outfilemaxmin.close();
	
		if(learning_method == "bma" && grid_scaleM == 5) {
			string datafile = ofilenamerep+"data.txt";
			// Out put simulated data for reproducibility purposes
			// Only do this when the method is bma
			ofstream outfiledata(datafile.c_str());
			for(int i=0; i < p.samplesize; i++) {
				for(int d=0; d<(int)p.data[0].size();d++) {
					outfiledata<<p.data[i][d]<<"\t";
				}
				outfiledata<<endl;
			}
			for(int i=0; i < p.valisize; i++) {
				for(int d=0; d<(int)p.validata[0].size();d++) {
					outfiledata<<p.validata[i][d]<<"\t";
				}
				outfiledata<<endl;
			}
			outfiledata.close();			
		}
		
		double KLdist_shift =0.0;
		double KLdist = 0.0;
		vector<OneDensity_data> ensemble;	
		
		if(learning_method == "bma" || learning_method == "gibbs") {
			// fit Bayesian density forests
			timer.Begin();
			
			ensemble.clear();
			bma_density_forest(p, para, ensemble, ensemble_sizeP, n_shiftM, shift_sizeM, n_shiftP, shift_sizeP, betaM, betaP, ofilenamerep);
		
			timer.End();
			outTime_BSP_shift<<timer.RunTimeProfile()<<endl;
	
			// calculating KL divergence, linxi
			KLdist_shift =0.0;
		
			KLdist_shift = simulation_KL_fromf_shift(p.testdata, ensemble, data_model);
		
			outKL_BSP_shift<<KLdist_shift<<endl;
			cout<<"KL divergence from the true is (with random shifting): "<<KLdist_shift<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			
			// Fit optional polya tree without shift, linxi
			timer.Begin();
			// partitions for marginal dist's, linxi
		    vector<OnePartition_data> marginalP; // structure defined in sampling.h, linxi
			marginalP.clear();
			// tree for the joint distribution
			OnePartition_data bestlevelPs;
		
			optional_polya_tree_copula(p, para, marginalP, bestlevelPs, betaM, betaP, ofilenamerep);
			timer.End();
			outTime_BSP<<timer.RunTimeProfile()<<endl;
		    //+ end by ed520
	
			// calculating KL divergence, linxi
			KLdist = 0.0;
			KLdist = simulation_KL_fromf(p.testdata, bestlevelPs, marginalP, data_model);
		
			outKL_BSP<<KLdist<<endl;
			cout<<"KL divergence from the true is (BMA with random shifting): "<<KLdist_shift<<endl;
			cout<<"KL divergence from the true is (OPT): "<<KLdist<<endl;	
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			
		}
		
		if(learning_method == "gibbs") {
			// For Gibbs density forests
			timer.Begin();
		
			int depth_max = (int) (log ((double)(p.samplesize+p.valisize))) /(log(2.0));
			int depth_min = 3; //(int) ((0.2*log((double)p.samplesize))/(log(2.0))+0.2+2);
			int number_of_forests = 3*(depth_max - depth_min+1);
		
			// For testing purposes 
			//int depth_max = 6;
			//int depth_min = 6; 
			//int number_of_forests = 1;
			
			vector<OneForest_data> forest_of_Fs;	
			forest_of_Fs.clear();	
		
			// fit gibbs density forests
			gibbs_density_forest(p, para, forest_of_Fs, depth_max, depth_min, ensemble_sizeP, number_of_forests, betaM, betaP, splitP,ofilenamerep, tree_temp, forest_temp, Cdepth, false);
		
			// Find the MAP estimator, linxi
			int bestFind =0;
			double bestFpost_p =0;
			for(int l=0; l<(int)forest_of_Fs.size(); l++) {
				//cout<<"Forest "<<l<<" Posterior: "<<forest_of_Fs[l].post_p<<" Wegith: "<<forest_of_Fs[l].weight<<" Tree weight: "<<forest_of_Fs[l].tree_weight<<endl;
				if(forest_of_Fs[l].post_p > bestFpost_p) {
					bestFpost_p = forest_of_Fs[l].post_p;
					bestFind =l;
				}
			}
			
			timer.End();
			outTime_BSP_fof<<timer.RunTimeProfile()<<endl;
		
			vector<double> KLdist_GDF = {0.0, 0.0};
		
			KLdist_GDF = simulation_KL_fromf_fof(p.testdata, forest_of_Fs, data_model, bestFind);
			double KLdist_Gibbs = KLdist_GDF[0]; // KL for posterior mode
			double KLdist_fof = KLdist_GDF[1];  // KL for posterior mean
			outKL_BSP_fof<<KLdist_fof<<endl;
			outKL_BSP_Gibbs<<KLdist_Gibbs<<endl;
			cout<<"KL divergence from the true is (BMA with random shifting): "<<KLdist_shift<<endl;
			cout<<"KL divergence from the true is (OPT): "<<KLdist<<endl;	
			cout<<"KL divergence from the true is (with random shifting, applying FOF): "<<KLdist_fof<<endl;
			cout<<"KL divergence from the true is (with random shifting, Gibbs posterior mode): "<<KLdist_Gibbs<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;		
			
			for(int l=0; l<(int)forest_of_Fs.size(); l++) {
				cout<<"Forest "<<l<<" Posterior: "<<forest_of_Fs[l].post_p<<" Weight: "<<forest_of_Fs[l].weight<<" Tree weight: "<<forest_of_Fs[l].tree_weight<<endl;
			}	
			
			if(UQ_flag) {
				// UQ for Gibbs density forests
				// Only do this once in the last repetition 
				timer.Begin();
		
				if(iter == nrep-1) {
					string ofilenameUQ = ofilename + "UQ/";
					string clineUQ = "mkdir -p "+ofilenameUQ;
					const int dirUQ= system(clineUQ.c_str());
					if (dirUQ< 0) return 5;
		
					int depth_max = (int) (log ((double)(p.samplesize+p.valisize))) /(log(2.0));
					int depth_min = 3; //(int) ((0.2*log((double)p.samplesize))/(log(2.0))+0.2+2);
			
					vector<OneForest_data> GDF_posterior;	
					GDF_posterior.clear();
		
					// fit gibbs density forests
					// Draw 2000 forest topologies from the posterior 
					gibbs_density_forest(p, para, GDF_posterior , depth_max, depth_min, ensemble_sizeP, 2000, betaM, betaP, splitP, ofilenameUQ, tree_temp, forest_temp,  Cdepth, UQ_flag);
		
					// Sample from the posterior 
					vector<OneForest_data> UQ_samples;
					UQ_samples.clear();
			
					// Sample 200 forests from the posterior
					GDF_UQ_sampling(GDF_posterior, UQ_samples, 200, 0.5);
			
					// construct test points
					int total_num = 1024*4;
					vector<vector<double> > test_points;
					test_points.clear();
					test_points.resize(total_num);
					for(int i=0; i < total_num; i++) {
						test_points[i].clear();
						test_points[i].resize(p.dim);
						for(int d=0; d< p.dim; d++) {
							test_points[i][d] = (double)i/ (double)(total_num) + (double)1/(double)(2*total_num);
						}
					}
			
					// evaluate densities at test points
					GDF_UQ_path_evaluation(UQ_samples, test_points, ofilenameUQ);		
				}
			
				timer.End();
				outTime_UQ<<timer.RunTimeProfile()<<endl;
			}
			
			// Checking the depth
			cout<<"Maximum depth is: "<<depth_max<<"; Minimum depth is: "<<depth_min<<endl;	
			
			/*
			int temp_tree_size = (int)forest_of_Fs[0].forestP[0].marginalP[0].sregs.size();
			for(int i=0; i< temp_tree_size-1; i++) {
				cout<< "Cut "<<i<<" :"<<endl;
				cout<<forest_of_Fs[0].forestP[0].marginalP[0].cut_index[i]<<endl;
				cout<<forest_of_Fs[0].forestP[0].marginalP[0].left_count[i]<<endl;
				cout<<forest_of_Fs[0].forestP[0].marginalP[0].right_count[i]<<endl;	
			}
		
			vector<double> mle_density = get_all_density(forest_of_Fs[0].forestP[0].marginalP[0]);
			vector<double> sample_density = sample_posterior_density(forest_of_Fs[0].forestP[0].marginalP[0], 0.5);
			for(int i=0; i<temp_tree_size; i++) {
				cout<<mle_density[i]<<"\t"<<sample_density[i]<<endl;
			}
			*/
		}
    					
		if(learning_method == "RF") {
			// fit density random forests
			timer.Begin();
			
			//vector<OneDensity_data> ensemble;	
			ensemble.clear();
			density_RF(p, para, ensemble, ensemble_sizeP, n_shiftM, shift_sizeM, n_shiftP, shift_sizeP, betaM, betaP,  ofilenamerep, RFdepth);
		
			timer.End();
			outTime_RF<<timer.RunTimeProfile()<<endl;
	
			// calculating KL divergence, linxi
			double KLdist_RF=0.0;
		
			KLdist_RF = simulation_KL_fromf_shift(p.testdata, ensemble, data_model);
		
			outKL_RF<<KLdist_RF<<endl;
			cout<<"KL divergence from the true is (BMA with random shifting): "<<KLdist_shift<<endl;
			cout<<"KL divergence from the true is (OPT): "<<KLdist<<endl;	
			//cout<<"KL divergence from the true is (with random shifting, applying FOF): "<<KLdist_fof<<endl;
			//cout<<"KL divergence from the true is (with random shifting, Gibbs posterior mode): "<<KLdist_Gibbs<<endl;
			cout<<"KL divergence from the true for RF is: "<<KLdist_RF<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
			cout<<"-------------------------------------------------------------------------------"<<endl;
		
		}
	
	}	
	
	outKL_BSP.close();
	outKL_BSP_shift.close();
	outKL_BSP_fof.close();
	outKL_BSP_Gibbs.close();
	outTime_BSP.close();
	outTime_BSP_shift.close();
	outTime_BSP_fof.close();
	outTime_UQ.close();	
    return 0;

}

void print_usage_and_exit() {
    //    cerr << "Usage: " + c::PROG_NAME + " <BSP>" << "\n";
    cerr << "Options:" << "\n";
    cerr << "-== Density Estimation ==-" << '\n';
    cerr << "  BSP_C       -- MAP partitions for marginal distribution and MAP partition for joint distribution after copula transformation" << "\n";
    cerr << "  BSP_NC      -- MAP partition without copula" << "\n";
    cerr << "\n";
    //    cerr << "-== Other tools ==-" << '\n';
    //    cerr << "  hell_dist   -- Compute sample Hellinger distance from a known density" << "\n";
    //    cerr << "  density     -- Get the density at particular points" << "\n";
    exit(2);
}



int BSP_tree_C(vector<string> params ){
    string wrong_message = "ED520 tree version : wrong input";
    if(params.size()> 4 || params.size() < 2)
    {
        cout<<wrong_message<<endl;
        return 3;
    }
    else{
        cerr << "BSP_Tree is programmed by Parallel Computing System Lab (ED520) in NCTU\n";
        cerr << "https://sites.google.com/a/g2.nctu.edu.tw/parallel-computing-lab/" << "\n";
        cerr << "\n";
    }
    tree BspTree;
    BspTree.discrete=true;
    string ofilename = params[1] ;
    vector<double> mmax, mmin;    
    BspTree.Data = read_data(params[0], true, mmax, mmin);
    BspTree.Dim = BspTree.Data[0].size();
    if(params.size()==4)
        BspTree.Levels= strTo<int>(params[2]);
    else 
        BspTree.Levels = 1000;  // default # of cut
    if(params.size()==4)
        BspTree.n = strTo<int>(params[3]); 
    else 
        BspTree.n = 200;       // default # of partition
    
    string tempmaxmin= ofilename+"_maxmin.txt";
    ofstream outfilemaxmin (tempmaxmin.c_str());
    for(int d=0; d<BspTree.Dim   ;d++){
        outfilemaxmin << mmax[d]<<'\t'<< mmin[d]<<'\n';
    }
    outfilemaxmin.close();
    BspTree.mmax = mmax;	//+ 412
    BspTree.mmin = mmin;	//+ 412   
    
    BspTree.tree_nonpara_init( );
    // copula initialization : 
    tree BspTree_para;
    BspTree_para.discrete=true;
    BspTree_para.maxpercentage=BspTree.maxpercentage; //some parameter for statistics
    BspTree_para.Levels=BspTree.Levels;     
    BspTree_para.n=BspTree.n;          
    BspTree_para.Resampling=BspTree.Resampling;
    BspTree_para.Samplesize=BspTree.Samplesize;   
    BspTree_para.Steps=BspTree.Steps;
    BspTree_para.beta=BspTree.beta;
    // end of copula parameter !!! 
    string temptrans = ofilename+"_trans.txt";
    ofstream outfiletrans (temptrans.c_str());
    // for final output file 
#ifdef OUTPUT_FILE
    vector<OnePartition_data> marginalP;
    vector<pair <double , double> > temp_intervals;
    vector <vector<pair <double , double> > > marginal_intervals;
#endif
    vector <double> onevariable(BspTree.Samplesize,0);
    BspTree_para.para_share_init();
    int reg_count[BspTree.Dim], reg_count2[BspTree.Dim];		//+ 412 for testing
    //start Copula
    for(int i=0;i<BspTree.Dim;i++){
        cout<<"dim is "<<i<<endl;
        BspTree_para.tree_para_init(  );
        
        for(int j=0 ; j < BspTree_para.Samplesize ; j++ ){
            BspTree_para.Data_array[j]= BspTree.Data_array[j*BspTree.Dim+i];
            onevariable[j]=BspTree.Data_array[j*BspTree.Dim+i]; //  Data[j][i]; 
        }
        for(int j =0 ; j <BspTree.Samplesize; ++j ){
            BspTree_para.Point_index[j]= BspTree_para.Dim*j;    //  create all index first 
        }

        BspTree_para.SIS_tree();
       
        BspTree_para.List_index=BspTree_para.List_index_best;
#ifdef DEBUG
        BspTree_para.print_partition1_t(  BspTree_para.P_bestlevel );
#endif
        vector<double> newx = BspTree_para.Ftransform_t(onevariable,  BspTree_para.P_bestlevel);

        for(int j=0; j<10; j++) cerr<<onevariable[j]<<" ";
        cerr<<endl;

        for(int j=0; j<10; j++) cerr<<newx[j]<<" ";
        cerr<<endl;
#ifdef DEBUG
        vector<double> recoveronevariable = BspTree_para.inv_Ftransform_t(newx, BspTree_para.P_bestlevel );

        for(int j=0; j<10; j++) cerr<<recoveronevariable[j]<<" ";
        cerr<<endl;
#endif
        for(int j=0; j<BspTree.Samplesize; j++){
            BspTree.Data_array[j*BspTree.Dim+i] = newx[j];  //+ by ed520  need to be replaced 
            outfiletrans << newx[j]<<'\t';
        }
        outfiletrans<<'\n';
#ifdef OUTPUT_FILE
        string temp= ofilename+'_'+toStr<int>(i)+".txt";
        ofstream outfile (temp.c_str());
        BspTree_para.print_partition2_t( outfile, BspTree_para.P_bestlevel );
        temp_intervals = BspTree_para.onedimpartitiontointerval_t_copy( BspTree_para.P_bestlevel );
        marginal_intervals.push_back(temp_intervals);
        outfile.close();
        //BspTree_para.SaveTree();			//+ 412: Resulting exception
        reg_count[i] = BspTree_para.P_bestlevel.region[BspTree_para.List_index].size(); //+ 412: Testing
        reg_count2[i]= BspTree_para.Bsp_Tree.size(); 	//+ 412: Testing
        temp = ofilename+'_'+toStr<int>(i)+"_Tree.txt";	//+ 412
        ofstream outfile_t (temp.c_str());  		//+ 412
        BspTree_para.output_tree(outfile_t,false);	//+ 412
        outfile_t.close();				//+ 412
#endif
        for(unsigned del=0; del < BspTree_para.Bsp_Tree.size(); del++ )
        {
            BspTree_para.Bsp_Tree[del]->destroy(); 
        }
        cout<<" Finish on Dim of copula Dim is  "<<i<<"  \n\n";
    }
    //end Copula
    BspTree_para.tree_para_destroy();
    outfiletrans.close();
    
    //start Non-copula
    BspTree.SIS_tree();
    BspTree.List_index=BspTree.List_index_best;
    BspTree.SaveTree();		//modified 412
    
    // for output file 
#ifdef OUTPUT_FILE
    string temp= ofilename+"_Big.txt";
    string temp2 = ofilename+"_OriBig.txt";	//+ 412, records c(u1,...,uD)
    ofstream outfile (temp.c_str());
    ofstream outfile2(temp2.c_str());		//+ 412
    Partition &op = BspTree.P_bestlevel;
    
    for (unsigned i = 0; i < op.region[BspTree.List_index].size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < BspTree.Dim; d++) {
            range = convert_ranges(op.region[BspTree.List_index][i]->Region_mask[d]);
            double ub,lb;
            ub = BspTree.inv_Ftransform_for_output(range.second, marginal_intervals[d]  );
            lb = BspTree.inv_Ftransform_for_output(range.first,  marginal_intervals[d]  );
            outfile << lb << " " << ub <<" ";
            outfile2<< range.first << " " << range.second <<" ";
        }
        outfile << setprecision(12) <<(double)op.region[BspTree.List_index][i]->length/(double)BspTree.Samplesize/ exp(lprod_usint_mask(op.region[BspTree.List_index][i]->Region_mask)) <<" "<< op.region[BspTree.List_index][i]->length<<endl;
        outfile2 << setprecision(12) <<(double)op.region[BspTree.List_index][i]->length/(double)BspTree.Samplesize/ exp(lprod_usint_mask(op.region[BspTree.List_index][i]->Region_mask)) <<" "<< op.region[BspTree.List_index][i]->length<<endl;
    }
    outfile.close();
    outfile2.close();
#endif

    // end of output file 
    // end of process
    // end of the showing !! 

    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<BspTree.Dim<<endl;
    cout<<"Data_num = "<<BspTree.Samplesize<<endl;
    cout<<"# of Partition = "<<BspTree.n<<endl;   
#ifdef PROFILE
    cout<<"====== Profile Info ======"<<endl;
    cout<<"NC CPU COUNT = "<<timer.CPU_COUNT<<endl;
    cout<<"NC HIT COUNT = "<<timer.HIT_COUNT<<endl;
    cout<<"Memory Usage = "<<BspTree.Last_point_index*4/1024/1024<<" MB"<<endl;
    cout<<"Total Memory = "<<(long)(BspTree.n*BspTree.Samplesize*SCALE)*4/1024/1024<<" MB"<<endl;
#endif
    //BspTree.SaveTree();	//modified 412
    temp = ofilename+"_Tree.txt";
    ofstream outfile3(temp.c_str());
    BspTree.output_tree(outfile3,false);
    outfile3.close(); 
    
    cout << reg_count[0] <<endl << reg_count[1] <<endl;
    cout << reg_count2[0] <<endl<< reg_count2[1]<<endl;
    
    return 0;

}



int BSP_tree_NC(vector<string> params ){
    string wrong_message = "ED520 tree version : wrong input";
    if(params.size()> 4 || params.size() < 2)
    {
        cout<<wrong_message<<endl;
        return 3;
    }
    else{
        cerr << "BSP_Tree is programmed by Parallel Computing System Lab (ED520) in NCTU\n";
        cerr << "https://sites.google.com/a/g2.nctu.edu.tw/parallel-computing-lab/" << "\n";
        cerr << "\n";
    }
    tree BspTree;
    BspTree.discrete=true;
    string ofilename = params[1] ;
    vector<double> mmax, mmin;
    BspTree.Data = read_data(params[0], true, mmax, mmin);
    BspTree.mmax = mmax;	//+ 412
    BspTree.mmin = mmin;	//+ 412    
    BspTree.Dim = BspTree.Data[0].size();
    if(params.size()==4)
        BspTree.Levels= strTo<int>(params[2]);
    else
        BspTree.Levels = 1000;  // default # of cut
    if(params.size()==4)
        BspTree.n = strTo<int>(params[3]);
    else
        BspTree.n = 200;       // default # of partition
    BspTree.tree_nonpara_init( );
        
    BspTree.SIS_tree();
    BspTree.List_index=BspTree.List_index_best;
    
    BspTree.SaveTree();        //modified 412
#ifdef OUTPUT_FILE
    string temp= ofilename+"_Big.txt";
    ofstream outfile (temp.c_str());
    BspTree.print_partition3_t(outfile, BspTree.P_bestlevel, mmax, mmin);

    outfile.close();
#endif

    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<BspTree.Dim<<endl;
    cout<<"Data_num = "<<BspTree.Samplesize<<endl;
    cout<<"# of Partition = "<<BspTree.n<<endl;
#ifdef PROFILE
    cout<<"====== Profile Info ======"<<endl;
    cout<<"NC CPU COUNT = "<<timer.CPU_COUNT<<endl;
    cout<<"NC HIT COUNT = "<<timer.HIT_COUNT<<endl;
    cout<<"Memory Usage = "<<BspTree.Last_point_index*4/1024/1024<<" MB"<<endl;
    cout<<"Total Memory = "<<(long)(BspTree.n*BspTree.Samplesize*SCALE)*4/1024/1024<<" MB"<<endl;
#endif
    //BspTree.SaveTree();   			//modified 412
    
    temp = ofilename+"_Tree.txt";		//+ 412
    
    ofstream outfile2(temp.c_str());		//+ 412
    BspTree.output_tree(outfile2,true);		//+ 412
    outfile2.close();				//+ 412
    
    //test_coding(BspTree, ofilename);
    
    return 0;
}



int CDF(vector<string> & params){
    if (params.size() != 4) {
        cerr<<params.size()<<endl;
        cerr << "usage_text" << endl;
        return 3;
    }

    vector<vector<double> > test_data = read_data(params[1]);  //data

    int test_N = (int)test_data.size();    
    int dim    = (int)test_data[0].size();
    vector<double> CDF(test_N, 0);

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    if (params[0]=="-c"){
        vector<vector<double> > Fx = test_data; // +412 size is test_N*D
        vector<vector<vector<double> > > MarginalPs;
        vector<vector<double> > mmax_min = read_data( params[2] + "_maxmin.txt");  // params[2] for the partition_info file
        
        for (int i=0; i<test_N; i++){
            for (int d=0; d<dim; d++){
                test_data[i][d] = (test_data[i][d] - mmax_min[d][1]) / (mmax_min[d][0] - mmax_min[d][1]);
            }
        }
        for(int d=0; d<dim; d++){
            MarginalPs.push_back(read_data( params[2]+'_'+toStr<int>(d)+".txt"));
        }

        for (int i = 0; i < test_N; i++){
            for (int d = 0; d < dim; d++){	// Transform from x to Fx(x) for each dimension
                vector<double> onepiece;
                onepiece.push_back(test_data[i][d]);                
                Fx[i][d] = CDF_from_partition(onepiece, MarginalPs[d]);
            }
        }                                                
        vector<vector<double> > P = read_data( params[2]+"_OriBig.txt");	// Using C(u1,...,uD) instead of C(x1,...,xD)
        for (int i = 0; i < test_N; i++) {
            CDF[i] = CDF_from_partition(Fx[i], P);  //the 2*dim col of P should be density
        }
    }
    else  if (params[0]=="-n"){
        vector<vector<double> > P = read_data( params[2]+"_Big.txt");// params[2] for the partition_info file
        for (int i = 0; i < test_N; i++) {
            CDF[i] = CDF_from_partition(test_data[i], P);
        }
    }
    else {
        cerr << "usage_text" << endl;
        return 1;
    }

    string tempout= params[3]+"_CDF.txt";   
    ofstream outfileout (tempout.c_str());
    for(int i=0; i<test_N; i++){
        outfileout <<setprecision(8)<<CDF[i]<<'\n';
    }

    outfileout.close();
        
    return 0;

}



int CDF_tree(vector<string> & params){
    if (params.size() != 4) {
        cerr<<params.size()<<endl;
        cerr << "usage_text" << endl;
        return 3;
    }

    vector<vector<double> > test_data = read_data(params[1]);  //data

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();
    vector<double> CDF(test_N, 0);

    cerr << test_N << " data points in " << dim << " dimensions.\n";

    if (params[0]=="-c"){
        vector<vector<double> > Fx = test_data; // +412 size is test_N*D
        vector<vector<vector<double> > > MarginalPs;
        vector<vector<double> > mmax_min = read_data( params[2] + "_maxmin.txt");  // params[2] for the partition_info file
        for (int i=0; i<test_N; i++){
            for (int d=0; d<dim; d++){
                test_data[i][d] = (test_data[i][d] - mmax_min[d][1]) / (mmax_min[d][0] - mmax_min[d][1]);
            }
        }
        for(int d=0; d<dim; d++){
            MarginalPs.push_back(read_data( params[2]+'_'+toStr<int>(d)+".txt"));
        }

        for (int i = 0; i < test_N; i++){
            for (int d = 0; d < dim; d++){				// Transform from x to Fx(x) for each dimension
                vector<double> onepiece;
                onepiece.push_back(test_data[i][d]);                
                Fx[i][d] = CDF_from_partition(onepiece, MarginalPs[d]);
            }
        }                                                
        //vector<vector<double> > P = read_data( params[2]+"_OriBig.txt");// Using C(u1,...,uD) instead of C(x1,...,xD)
        vector<vector<double> > P = read_data( params[2]+"_Tree.txt");// Using C(u1,...,uD) instead of C(x1,...,xD)
        for (int i = 0; i < test_N; i++) {
            //CDF[i] = CDF_from_partition(Fx[i], P);
            CDF[i] = CDF_from_tree(Fx[i], P);
            //cout<< "Fu(u1,...,uD) = " << CDF[i] <<" (From Partition)"<<endl;
            //cout<< "Fu(u1,...,uD) = " << CDF_from_tree(Fx[i], P) <<" (From Tree)"<<endl;
        }
    }
    else  if (params[0]=="-n"){
        vector<vector<double> > P = read_data( params[2]+"_Tree.txt");	// params[2] for the partition_info file
        cout<<endl;    
        for (int i = 0; i < test_N; i++) {
            CDF[i] = CDF_from_tree(test_data[i], P);
        }
    }
    else {
        cerr << "usage_text" << endl;
        return 1;
    }

    string tempout= params[3]+"_CDF.txt";   
    ofstream outfileout (tempout.c_str());
    for(int i=0; i<test_N; i++){
        outfileout <<setprecision(8)<<CDF[i]<<'\n';
    }

    outfileout.close();

    return 0;
}



void test_partition_generation( Partition &mytestpartition ,  Region_Node *  mytestregion_a    ){

    int paradim=1;
    for( int ii=0 ; ii <10 ; ii ++){
        mytestregion_a[ii].length = ii;
        //mytestregion_b[i].length = 0;
        if(ii==0){
            mytestregion_a[ii].um.mask=0;
            mytestregion_a[ii].um.x=0;
            for(int j= 0 ; j <paradim ; j++)
                mytestregion_a[ii].Region_mask.push_back(mytestregion_a[ii].um);
        }
        else{
            for(int j = 0 ; j < paradim ; j ++){
                mytestregion_a[ii].Region_mask=mytestregion_a[ii-1].Region_mask;
                mytestregion_a[ii].Region_mask[j].x = mytestregion_a[ii-1].Region_mask[j].x<<1;
                mytestregion_a[ii].Region_mask[j].mask = (mytestregion_a[ii-1].Region_mask[j].mask<<1)+1 ;
            }
        }
    }
    for(int ii=0 ; ii <10 ; ii++)
    {
        mytestpartition.region_1.push_back(& mytestregion_a[ii]);
    }
    mytestpartition.region_2.push_back(NULL);
    mytestpartition.region.push_back(mytestpartition.region_1);
    mytestpartition.region.push_back(mytestpartition.region_2);


}

void test_coding( tree	BspTree, string ofilename ){
    string temp = ofilename+"_Partition.txt";
    ofstream outfile (temp.c_str());
    BspTree.print_partition2_t(outfile,BspTree.P_bestlevel);
    outfile.close();
 
    cout<<"#Data = "<<BspTree.Samplesize<<endl;
    //BspTree.P_bestlevel.Samplesize = BspTree.Samplesize;
    //cout<< BspTree.P_bestlevel.Samplesize << endl;

    double Prob, sumProb = 0;
    for (unsigned int i=0; i<BspTree.P_bestlevel.region[BspTree.List_index].size(); i++){
        Prob = BspTree.P_bestlevel.region[BspTree.List_index][i]->prob;
        sumProb += BspTree.P_bestlevel.region[BspTree.List_index][i]->prob;        
    }
    cout<<"Total Prob: "<< sumProb <<endl;                        
    
    pair<double,double> reg_range = convert_ranges(BspTree.P_bestlevel.region[BspTree.List_index][0]->um);
    
    cout<<"Size of Region_mask (dim): "<<BspTree.P_bestlevel.region[BspTree.List_index][0]->Region_mask.size()<<endl;    
    for (int i=0; i<BspTree.Dim; i++){
        reg_range = convert_ranges(BspTree.P_bestlevel.region[BspTree.List_index][0]->Region_mask[i]);
        cout<<"Lower bound: "<<reg_range.first<<" "<<"Upper bound: "<<reg_range.second<<endl;
    }   
    
    for (unsigned int i=0; i<BspTree.P_bestlevel.region[BspTree.List_index].size();i++){
        //cout << "Volume: " << BspTree.P_bestlevel.region[BspTree.List_index][i]->volume <<'\t';
        //cout << "Density: " << BspTree.P_bestlevel.region[BspTree.List_index][i]->prob/BspTree.P_bestlevel.region[BspTree.List_index][i]->volume <<'\t';
        //cout << "#data: " << BspTree.P_bestlevel.region[BspTree.List_index][i]->length << endl;
    }
    //
    cout<<"============================================"<<endl;
    Region_Node* node = BspTree.P_bestlevel.region[BspTree.List_index][15];
    cout<<"Current node's #data: "<<'\t'<< node->length << endl;
    cout<<"Its child size: "<< node->Node_child.size() << endl;
    
    for (int d=0; d<BspTree.Dim; d++){
        reg_range = convert_ranges(node->Region_mask[d]);
        cout<< reg_range.first << '\t' << reg_range.second <<'\t';
    }
    cout<<endl;    
    
    cout<<"Parent node's #data: "<<'\t'<< node->Node_parent->length << endl;
    cout<<"Its child size: "<< node->Node_parent->Node_child.size() << endl;
    cout<<"------->>>>> Cut dimension: "<< node->Node_parent->cut_dim() <<" ------"<< endl;

    for (int d=0; d<BspTree.Dim; d++){
        reg_range = convert_ranges(node->Node_parent->Region_mask[d]);
        cout<< reg_range.first << '\t' << reg_range.second <<'\t';
    }
    cout<<endl;

    cout<<"Brother node's #data: "<<'\t'<< node->Node_parent->Node_child[1]->length << endl;
    cout<<"Its child size: "<< node->Node_parent->Node_child[1]->Node_child.size() << endl;

    for (int d=0; d<BspTree.Dim; d++){
        reg_range = convert_ranges(node->Node_parent->Node_child[1]->Region_mask[d]);
        cout<< reg_range.first << '\t' << reg_range.second << '\t';
    }
    cout<<endl;
    
    cout<<"node->Node_parent->Node_child[0]'s #data: "<< node->Node_parent->Node_child[0]->length<<endl;
    cout<<"node->Node_parent->Node_child[1]'s #data: "<< node->Node_parent->Node_child[1]->length<<endl;
    
    while (true){
        double amount;
        node = node->Node_parent;
        amount = node->length;
        cout<<"#Data = "<< amount <<endl;
        if (amount == BspTree.Samplesize){
            break;
        }
    }
        
    cout<<'\t'<<" BSP Tree's root: "<< "#data = "<<BspTree.Root.length<<" "<<endl;
    for (int d=0; d<BspTree.Dim; d++){
        reg_range = convert_ranges(BspTree.Root.Region_mask[d]);
        cout<< " " << reg_range.first << '\t' << reg_range.second << '\t';
    }  
    cout<<endl;
    
    // Test Display
    temp = ofilename+"_Tree.txt";
    vector < vector <double> > Tree_Info = read_data(temp);
    
    cout<<"=========================================================================================="<<endl;
    cout<<"ID "<<'\t';
    for (int i=0; i<BspTree.Dim; i++){
        cout<<"LB"<<'\t'<<"UB"<<"\t";
    }
    cout<<"Density"<<'\t'<<"#Data"<<'\t'<<"CutDim"<<'\t'<<"Parent"<<'\t'<<"Child_1"<<'\t'<<"Child_2"<<endl;
    cout<<"------------------------------------------------------------------------------------------"<<endl;
    
    for (unsigned int i=0; i<20; i++){
        for (unsigned int j=0; j<Tree_Info[0].size(); j++){
            cout<< setiosflags(ios::right)<< Tree_Info[i][j] <<'\t';
        }
        cout<<endl;
    }
    //
    Region_Node* Node = BspTree.P_bestlevel.region[BspTree.List_index][0];
    pair<double,double> ranges;
    usint_mask	t_mask;
    t_mask.x = 0;
    t_mask.mask = 0;
    ranges = convert_ranges(t_mask);
    cout<< "Convert_ranges(all 0): "<<  ranges.first <<'\t'<< ranges.second <<endl;
    for (int d = 0; d < BspTree.Dim; d++){
        cout<< "Dim " << d <<": "<<endl;
        cout<< "Region_mask.x (Decimal): " << "Region_mask.mask (Decimal): " <<endl;
        cout<< Node->Region_mask[d].x <<'\t'<< Node->Region_mask[d].mask <<endl;
        cout<< dec2bin(Node->Region_mask[d].x) <<'\t'<< dec2bin(Node->Region_mask[d].mask) <<endl;
        ranges = convert_ranges(Node->Region_mask[d]);
        cout<< "Real ranges: "<<endl;
        cout<< ranges.first <<'\t'<< ranges.second <<endl;
    }
}
int dec2bin(int n){
    int binstream = 0;
    while (n>0){
        binstream *= 10;
        binstream += n - (n/2)*2;
        n = n/2;
    }
    return binstream;
}
