
#ifndef LETTER_H
#define	LETTER_H

#include "general_utils.h"
#include "output.h"
#include "SISfunctions.h"
#include "buildpathcttable.h"
#include "sampling.h"
#include "data_store.h"

// end_col_normalization=FALSE for classification problem
inline vector<vector<double> > read_data(string filename, bool end_col_normalization, vector<double>& mmax, vector<double>& mmin) {
    int dim = 0;
   cerr<<"start read data!"<<endl;
    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();

    infile.close();

    infile.open(filename.c_str());

    vector<vector<double> > data;

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
        data.push_back(d);
    }

    int k= dim;
    if(!end_col_normalization){ k=k-1;}

    vector<double> marginalmax(k, -100000000);
    vector<double> marginalmin(k,  1000000000);
    for(int i=0; i<k; i++){
        for (int j=0; j<(int)data.size(); j++){
            if(data[j][i]>marginalmax[i])  marginalmax[i]=data[j][i];
            if(data[j][i]<marginalmin[i])  marginalmin[i]=data[j][i];

        }
        cerr<<"max="<<marginalmax[i]<<"; min="<<marginalmin[i]<<endl;
    }
    for (int i = 0; i < k; i++) {
        if (marginalmin[i] == marginalmax[i]) {
            cout << "only one value in dim=" << i << endl;
            for (int j = 0; j < (int) data.size(); j++) {
                 data[j][i] = 0.45;
            }
            continue;
        }
        for (int j = 0; j < (int) data.size(); j++) {
            data[j][i] = (data[j][i] - marginalmin[i]) / (marginalmax[i] - marginalmin[i]);

        }
    }
	
	mmax = marginalmax;
	mmin = marginalmin;

    return data;
}


// end_col_normalization=FALSE for classification problem
inline vector<vector<double> > read_data(string filename) {
    int dim = 0;
   cerr<<"start read data!"<<endl;
    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();
    cerr<<"dim="<<dim<<endl;
    infile.close();

    infile.open(filename.c_str());

	// This is the data matrix
	vector<vector<double> > data;

    // Use the while loop to read in data
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
   //     cerr<<d[0]<<'\t';
        data.push_back(d);
    }

    return data;

}


// scale data to [a, b], also get the range of data
// end_col_normalization=FALSE for classification problem
// para "scale" indicates whether or not to scale data
inline void scale_data(vector<vector<double> >& data, double a, double b, bool end_col_normalization, bool scale, vector<double>& mmax, vector<double>& mmin) {
	int dim=0;
	dim = (int)data[0].size();
	
	int k= dim;
	if(!end_col_normalization){ k=k-1;}
	
    vector<double> marginalmax(k, -100000000);
    vector<double> marginalmin(k,  1000000000);
    for(int i=0; i<k; i++){
        for (int j=0; j<(int)data.size(); j++){
            if(data[j][i]>marginalmax[i])  marginalmax[i]=data[j][i];
            if(data[j][i]<marginalmin[i])  marginalmin[i]=data[j][i];

        }
        cerr<<"max="<<marginalmax[i]<<"; min="<<marginalmin[i]<<endl;
    }
	
	if(scale) {
	    int sumsmall=0;
	    for (int i = 0; i < k; i++) {
	        if (marginalmin[i] == marginalmax[i]) {
	            cout << "only one value in dim=" << i << endl;
	            for (int j = 0; j < (int) data.size(); j++) {
	                 data[j][i] = 0.45;
	            }    
	            continue;
	        }
	        for (int j = 0; j < (int) data.size(); j++) {
	  //          data[j][i] = (data[j][i] - marginalmin[i] + 0.001) / (marginalmax[i] - marginalmin[i] + 0.002);
	            data[j][i] = a + (b-a)*(data[j][i] - marginalmin[i]) / (marginalmax[i] - marginalmin[i]);
	            if(data[j][i]< a ) data[j][i]=a;
	            if(data[j][i]> b) data[j][i]= b;

	//            cerr<<data[j][i]<<'\t';
	            if(data[j][i]<a+0.01) sumsmall+=1;

	        }
	    }
	    cerr<<"sumsmall="<<sumsmall<<endl;	
	}
    
    mmax = marginalmax;
    mmin = marginalmin;
}

// scale data, multiply by 1/2 (downsize = T) or 2 (downsize =F), linxi
inline void scale_by_two(double* data, int samplesize, int dim, bool downsize) {
	double a;
	if(downsize) a = 0.5;
	else a = 2.0;
	
	for(int i=0; i < samplesize; i++) {
		for(int d=0; d< dim; d++) {
			data[i*dim+d] = a*data[i*dim+d];
		}
	}
}

// scale data, multiple by 1/2 or 2 input is a matrix, linxi 
inline void scale_by_two(vector<vector<double> >& data, bool downsize) {
	int samplesize = (int)data.size();
	int dim = (int)data[0].size();
	double a;
	if(downsize) a = 0.5;
	else a = 2.0;
	
	for(int i=0; i < samplesize; i++) {
		for(int d=0; d< dim; d++) {
			data[i][d] = a*data[i][d];
		}
	}
}

// shift data1D by a vector v, linxi
inline void shift_data(double* data, int samplesize, int dim, vector<double> shift) {
	for(int i =0; i<samplesize; i++) {
		for(int d=0; d<dim; d++) {
			data[i*dim+d] = data[i*dim+d] +shift[d];
		}
	}
}

// inverse transformation of shifting data by a vector v, linxi
inline void inv_shift_data(double* data, int samplesize, int dim, vector<double> shift) {
	for(int i =0; i<samplesize; i++) {
		for(int d=0; d<dim; d++) {
			data[i*dim+d] = data[i*dim+d] -shift[d];
		}
	}
}

inline vector<vector<double> > read_partition(string filename, bool end_line) {
    int dim = 0;

    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();

    infile.close();

    infile.open(filename.c_str());

    vector<vector<double> > data;

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

        bool good_data = true;
        for (int i = 0; i < dim; i++) {
            d[i] = strTo<double>(ll[i]);
            if (end_line) {
                if (i != (dim - 1) && (d[i] > 1.0 || d[i] < 0)) {   //parentheses by thchiu
                    good_data = false;
                    break;
                }
            } else {
                if (d[i] > (1.0 - (1.0 / 2147483648)) || d[i] < (0 + (1.0 / 2147483648))) {
                    good_data = false;
                    break;
                }
            }
        }
        if(good_data){
            data.push_back(d);
        }else{
            cerr << "Warning: Data out of range("<< dim <<"): " << line << '\n';
        }

    }
    return data;
}


#endif







