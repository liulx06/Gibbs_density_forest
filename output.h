/*
 * File:   output.h
 * Author: Mouse
 *
 * Created on November 2, 2011, 8:20 PM
 */

#ifndef OUTPUT_H
#define	OUTPUT_H

#include "sampling.h"

 inline void print_matrix( vector<vector<double> > & M){
     cout<<"distance matrix:\n";
     int m1 = (int)M.size();
     int m2 = (int)M[0].size();
     for(int i=0; i< m1; i++){
         for(int j=0; j<m2; j++){
             cout<<M[i][j]<<"\t";
         }
         cout<<"\n";
     }
 }


inline void print_partition(const OnePartition &P) {
    vector<shrink_reg> sr = P.sregs;
    for (int i = 0; i < (int) sr.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < (int) sr[0].reg_code.size(); d++) {
            range = convert_ranges(sr[i].reg_code[d]);
            cout << range.first << " " << range.second << " ";
            cout<<sr[i].reg_code[d].mask;
        }
        cout << 1 << endl;
    }
}

inline void print_partition(const OnePartition_data &P) {
//    vector<shrink_reg_data_short> sr = P.sregs;
    int totalnum = 0;

    for (int i = 0; i < (int) P.sregs.size(); i++) {
        totalnum += P.sregs[i].num;
    }


    for (int i = 0; i < (int) P.sregs.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < (int) P.sregs[0].reg_code.size(); d++) {
            range = convert_ranges(P.sregs[i].reg_code[d]);
            cout << range.first << " " << range.second << " ";
     //       cout<<sr[i].reg_code[d].mask<<" ";
        }
        cout <<  (double) P.sregs[i].num / totalnum / exp(lprod_usint_mask(P.sregs[i].reg_code)) << " " << P.sregs[i].num << endl;
    }
}


inline void print_partition(const OnePartition_data &P, const vector<double> & mmax, const vector<double>& mmin) {
//    vector<shrink_reg_data_short> sr = P.sregs;
    int totalnum = 0;

    for (int i = 0; i < (int) P.sregs.size(); i++) {
        totalnum += P.sregs[i].num;
    }


    for (int i = 0; i < (int) P.sregs.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < (int) P.sregs[0].reg_code.size(); d++) {
            range = convert_ranges(P.sregs[i].reg_code[d]);
            cout << range.first *(mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001 << " " << range.second *(mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001 << " ";
     //       cout<<sr[i].reg_code[d].mask<<" ";
        }
        cout <<  P.sregs[i].num<< (double) P.sregs[i].num / totalnum / exp(lprod_usint_mask(P.sregs[i].reg_code)) << " " << P.sregs[i].num << endl;
    }
}

inline void print_partition( ofstream& outfile, const OnePartition_data &P) {
    vector<shrink_reg_data_short> sr = P.sregs;
    int totalnum = 0;

    for (int i = 0; i < (int) sr.size(); i++) {
        totalnum += sr[i].num;
    }


    for (int i = 0; i < (int) sr.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < (int) sr[0].reg_code.size(); d++) {
            range = convert_ranges(sr[i].reg_code[d]);
             outfile << range.first << " " << range.second << " ";
        }
        outfile << (double) sr[i].num / totalnum / exp(lprod_usint_mask(sr[i].reg_code))<<" "<<P.sregs[i].num<<endl;

    }
}

inline void print_partition( ofstream& outfile, const OnePartition_data &P, const vector<double> & mmax, const vector<double>& mmin) {
    vector<shrink_reg_data_short> sr = P.sregs;
    int totalnum = 0;

    for (int i = 0; i < (int) sr.size(); i++) {
        totalnum += sr[i].num;
    }


    for (int i = 0; i < (int) sr.size(); i++) {
        pair<double, double> range;
        double newarea=1.0;
        for (int d = 0; d < (int) sr[0].reg_code.size(); d++) {
            range = convert_ranges(sr[i].reg_code[d]);
            double newrange1 = range.first *(mmax[d] - mmin[d])+ mmin[d];
            double newrange2 = range.second *(mmax[d] - mmin[d])+ mmin[d];
            newarea= newarea*(newrange2-newrange1);
             outfile <<newrange1 << " " <<newrange2 << " ";
        }
        outfile <<P.sregs[i].num/(double)totalnum/newarea<<" " << P.sregs[i].num<<endl;

    }
}




/*
inline double Fprob(double x, vector<pair<double,double> >& intervals){
    double y=0;
    bool cont=true;
    for(int i=0; cont && i<(int)intervals.size()-1; i++){
        if( x>intervals[i+1].first){
            y += (intervals[i+1].first-intervals[i].first)* intervals[i].second;
        }
        else {
            y += (x-intervals[i].first)* intervals[i].second;
            cont=false;
        }
    }
    if(cont) {
        int i=(int)intervals.size()-1;
        y += (x-intervals[i].first)* intervals[i].second;

    }
    if(y==0)  y=0.000001;                // make sure 0 and 1 don't appear.
    if(y>0.999999)  y=0.999999;
    return y;

}
*/

inline double Fprob(double x, vector<pair<double,double> >& intervals){
    double y=0;
    bool cont=true;
    if((int)intervals.size()==1) return x;
    if( x<intervals[1].first)  {
        y=intervals[0].second*x/intervals[1].first;
        cont=false;
    }
    for(int i=1; cont && i<((int)intervals.size()-1); i++){
        if( x<intervals[i+1].first){
            y = intervals[i-1].second+(intervals[i].second-intervals[i-1].second)*(x-intervals[i].first)/(intervals[i+1].first-intervals[i].first);
            cont=false;
        }
    }
    if(cont) {
        int i=(int)intervals.size()-1;
        if(1-intervals[i].first>1e-5){
                y = intervals[i-1].second+(1-intervals[i-1].second)*(x-intervals[i].first)/(1-intervals[i].first);
        }else y=intervals[i].first;
    }
    if(y>=-1e-5 && y< 1e-5)  y=1e-5;         // make sure 0 and 1 don't appear.
    else if(y>=1-1e-5 && y <= 1+1e-5)  y=1-1e-5;
    else if (y<-0.00001 || y>1.00001){
        cerr << "y: " << y << '\n';
        cerr << "wrong y in Fprob!";
        exit(10);
    }
    if(isnan(y)){
        cerr<<"y="<<y<<endl;
        //for(int i=0; i<intervals.size(); i++){    //- by ed520
        for(unsigned i=0; i<intervals.size(); i++){      //+ by ed520
            cerr<<intervals[i].first<<"\t"<<intervals[i].second<<endl;
        }
    }
    return y;

}


inline double inv_Fprob(double Fx, vector<pair<double,double> >& intervals){
    double x=-100000;
    if((int)intervals.size()==1) return Fx;
    if( Fx<intervals[0].second)  {
        x=intervals[1].first*Fx/intervals[0].second;
        return x;
    }
    for(int i=1; i<((int)intervals.size()-1); i++){
        if( Fx<intervals[i].second){
			if(intervals[i].second-intervals[i-1].second >1e-5) {
	            x = intervals[i].first+(intervals[i+1].first-intervals[i].first)*(Fx-intervals[i-1].second)/(intervals[i].second-intervals[i-1].second);
			}else {
				x = intervals[i].first;
			}
            return x;
        }
    }
    int i= (int)intervals.size()-1;
    if(1-intervals[i-1].second>1e-5){
        x = intervals[i].first+(1-intervals[i].first)*(Fx-intervals[i-1].second)/(1-intervals[i-1].second);
    }
    else x = intervals[i].first;
    return x;

}



inline bool compareinterval (const pair<double, double> &i, const pair<double, double> &j) {
	return (i.first<j.first);
}
// transform the partition to the pair<double,double> recording the starting point and the cdf. order the intervals by the position. the first pair.first=0
inline vector<pair<double,double> > onedimpartitiontointerval(OnePartition_data &P){
//    vector<shrink_reg_data_short> regs = P.sregs;
    int numreg = (int)P.sregs.size();
    pair<double,double> range;
    int totalnum=0;
    for(int i=0; i< numreg; i++){   totalnum += P.sregs[i].num;  }
    vector<pair<double,double> > intervals;
     if(numreg==1) {
        range.first=0;
        range.second=1;
        intervals.push_back(range);
        return (intervals);
    }
    for(int i=0; i< numreg; i++){
	range.first = convert_ranges(P.sregs[i].reg_code[0]).first;
        range.second = (double)P.sregs[i].num/totalnum/exp(lprod_usint_mask(P.sregs[i].reg_code)); //
//      range.second = intervals[i-1].second + (range.first- intervals[i-1].first)* range.second;
        intervals.push_back(range);
    }
    sort(intervals.begin(), intervals.end(), compareinterval);

   //  transform the intervals[i].second to the cdf.
    intervals[0].second =  (intervals[1].first- intervals[0].first)* intervals[0].second;
    for(int i=1; i<(numreg-1); i++){
        intervals[i].second = intervals[i-1].second + (intervals[i+1].first- intervals[i].first)* intervals[i].second;
    }
    intervals[numreg-1].second = intervals[numreg-2].second + (1- intervals[numreg-1].first)* intervals[numreg-1].second;

    if(fabs(intervals[numreg-1].second -1) >0.00001){
        cout<<"wrong in onedimpartitiontointerval\n";
        for(int i=0; i<numreg; i++){
            cout<<intervals[i].first<<" "<<intervals[i].second<<"\n";
        }
        exit(9);
    }

    return intervals;
}

 // return the indexes of the intervals overlap with range
 // notice: normalize ahead
 inline vector<pair<double,double> > search_overlap(pair<double, double> &range, vector<pair<double,double> >& interval){

     vector<pair<double,double> > trueintervals;
     pair<double, double> temp = range;
     if((int)interval.size()==1){
          trueintervals.push_back(temp);
          return trueintervals;
     }
     for(int i=0; i<(int)interval.size(); i++){
         if(interval[i].first>range.second){
             break;
         }
         if(interval[i].first>range.first){
             temp.second = interval[i].first;
             trueintervals.push_back(temp);
             temp.first= temp.second;
         }
     }
     temp.second = range.second;
     trueintervals.push_back(temp);
     return trueintervals;

 }

inline void output_full_partition( vector< vector<pair<double,double> > >& trueintervals, vector<vector<double> >& combinetrueintervals){
    combinetrueintervals.clear();

    vector<double> onetemp0(2,0);
    if((int)trueintervals.size() ==1){
         for(int i=0; i< (int)trueintervals[0].size(); i++){
             onetemp0[0] = trueintervals[0][i].first;
             onetemp0[1] = trueintervals[0][i].second;
             combinetrueintervals.push_back(onetemp0);
         }
        return;
    }

     vector<vector<double> > temp;
     vector<pair<double,double> > lastdim =  trueintervals[(int) trueintervals.size()-1];
    trueintervals.pop_back();
    output_full_partition( trueintervals,temp);
    for(int i=0; i< (int)lastdim.size(); i++){
        for(int j=0; j<(int)temp.size(); j++){
            temp[j].push_back(lastdim[i].first);
            temp[j].push_back(lastdim[i].second);
            combinetrueintervals.push_back(temp[j]);
            temp[j].pop_back();
            temp[j].pop_back();
        }
    }
    return;
}


// from the joint distribution partition and the one dim partitions, reconstruct the full partition.
 inline void recover_fullregions(  vector<pair<double, double> >& ranges,  vector<OnePartition_data> &marginalP,  vector<double>& mmax,vector<double> & mmin, ofstream  &outfile){
     int dim= (int)ranges.size();
     vector< vector<pair<double,double> > > intervals,trueintervals;
     vector<vector<double> > combined_trueintervals;
     // for each dimension, the intervals involved
     for(int d=0; d<dim; d++){
        intervals.push_back( onedimpartitiontointerval(marginalP[d]));
     }
     for(int d=0; d<dim; d++){
          trueintervals.push_back(search_overlap(ranges[d], intervals[d]));
 //         cerr<<"overlap for dim "<<d<<"="<<trueintervals[d].size()<<endl;
 //         cerr<<intervals[d].size()<<endl;
 //         cerr<<"ranges="<<ranges[d].first<<" "<<ranges[d].second<<endl;
     }
     output_full_partition( trueintervals, combined_trueintervals);
     cerr<<"combinedregion num="<<combined_trueintervals.size()<<endl;
//     double sumarea=0;
     for(int i=0; i<(int)combined_trueintervals.size(); i++){
  //       double onearea=1;

  //       for(int j=0; j<2; j++){
  //           onearea= onearea*(combined_trueintervals[i][2*j+1]-combined_trueintervals[i][2*j]);
  //       }
  //       sumarea=sumarea+onearea;

  //       cerr<<"sum="<<sumarea<<endl;

         for(int d=0; d<(int)combined_trueintervals[i].size()/2; d++){

             outfile << combined_trueintervals[i][2*d]*(mmax[d] - mmin[d])+ mmin[d]<<" ";
             outfile << combined_trueintervals[i][2*d+1]*(mmax[d] - mmin[d])+ mmin[d]<<" ";
         }
         outfile<<'\n';
     }
  //   cerr<<"areasum="<<sumarea<<endl;
     return;
 }

inline double Ftransform( double x, OnePartition_data &P){
    vector<pair<double,double> > intervals = onedimpartitiontointerval(P);

    return  Fprob(x, intervals);
}


// F* transformation for one dimensional data. return F(x) based on the OnePartition P
inline vector<double> Ftransform(vector<double> x, OnePartition_data &P){
    vector<double> F((int)x.size(), 0);
    vector<pair<double,double> > intervals = onedimpartitiontointerval(P);
    for(int i=0; i<(int)x.size(); i++)
    {
        if(i < 10){
            cerr << "Intervals: " << intervals[0].first
                    << "," << intervals[0].second  <<  '\n';
        }
        F[i] = Fprob(x[i], intervals);
    }
    return F;
}


inline void Ftransform( vector<vector<double> >& testdata, vector<vector<double> >& dataF, vector<OnePartition_data> &Ps){
    unsigned samplesize = testdata.size();
    unsigned dim = testdata[0].size();
    if(Ps.size()!= dim) {
        cout << "wrong number of OnePartition_data in Ftransform!\n";
        exit(8);
    }
    dataF = testdata;
    vector<double> onevariable(samplesize, 0);

    for (unsigned d = 0; d < dim; d++) {
        for (unsigned i = 0; i < samplesize; i++) {
            onevariable[i] = testdata[i][d];
        }
        vector<double> newx = Ftransform(onevariable, Ps[d]);
        for (unsigned i = 0; i < samplesize; i++) dataF[i][d] = newx[i];
    }
    return;
}


// F*^-1 transformation for one dimensional data. return F^-1(x) based on the OnePartition P
inline vector<double> inv_Ftransform(vector<double> Fx, OnePartition_data &P){
    vector<double> x((int)Fx.size(), 0);
    vector<pair<double,double> > intervals = onedimpartitiontointerval(P);
//    for(int i=0; i<(int)intervals.size(); i++) cerr<<intervals[i].first<<" "<<intervals[i].second<<endl;
    for(int i=0; i<(int)Fx.size(); i++)  x[i] = inv_Fprob(Fx[i], intervals);
    return x;
}

// F*^-1 transformation for one dimensional data. return F^-1(x) based on the OnePartition P
inline double inv_Ftransform(double Fx, OnePartition_data &P){

    vector<pair<double,double> > intervals = onedimpartitiontointerval(P);
//    for(int i=0; i<(int)intervals.size(); i++) cerr<<intervals[i].first<<" "<<intervals[i].second<<endl;
    double x= inv_Fprob(Fx, intervals);
    return x;
}

inline void write_sample(vector<vector<double> >data, string filename) {
    ofstream fs(filename.c_str());

    for (int i = 0; i < (int) data.size(); i++) {
        for (int j = 0; j < (int) data[i].size(); j++) {

            fs << data[i][j] << " ";
        }
        fs<<"\n";
    }
    fs.close();
}


// read data from file. If by row: The context of the file should be of the form: x1 y1,z1, x2,y2,z2,x3,y3,z3
// df=2 for just x and y. df=3 for x,y,z. samplesize=n.

inline void readdata(string fileName, int samplesize, int df, vector<vector<double> >& data, bool normalization) {
    cout << "samplesize = " << samplesize << endl;
    ifstream infile;

    infile.open(fileName.c_str());
    if (infile.fail()) {
        cout << "File could not be found." << endl;
        return;
    }

    vector<double> a(df, 0);
    vector<double> high(df, 0);
    vector<double> low(df, 0);
//    double r;
    for (int j = 0; j < samplesize; j++) {
        int i = 0;

        while (i < df) {
            infile >> a[i];
            if (a[i] > high[i]) high[i] = a[i];
            if (a[i] < low[i]) low[i] = a[i];
            i++;
        }
        if (j == 0) {
            low = a;
            high = a;
        }
 //       infile >> r; // ***********get rid of the extra variable
        data.push_back(a);
    }

     for (int i = 0; i < df; i++)   cout<<"high["<<i<<"]="<<high[i]<<"   low["<<i<<"]="<<low[i]<<endl;
    infile.close();
    if (normalization) {
        cout<<"start normalization\n";
        for (int j = 0; j < samplesize; j++) {
            for (int i = 0; i < df; i++) {

                data[j][i] = (data[j][i] - low[i]) / (high[i] - low[i]);
            }
        }
    }

/*    for (int j = 0; j < df; j++) {
        cout << j << "th dimension" << endl;
        for (int i = 0; i < samplesize; i++) {
            cout << data[i][j] << "\t";
        }
        cout << endl;
    }
*/
}


#endif	/* OUTPUT_H */

