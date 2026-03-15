


#ifndef SAMPLE_DENSITY_GENERATION_H
#define SAMPLE_DENSITY_GENERATION_H

#include <cmath>
#include "sampling.h"
#include "beta.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues" 


static const double L_PI=3.141596;

inline void sample_quarter_circle_2(int n, vector<vector<double> >& data, int d=2) {
	int dim = d;
	data.resize(n);
	int i = 0;
	while (i < n) {
		vector<double> temp;
		temp.resize(dim);
		for (int j = 0; j < dim; j++) {
			temp[j] = rand_double();
		}
		if (temp[d-1] * temp[d-1] + temp[d-2] * temp[d-2] > 1) {
			data[i] = temp;
			i++;
		} else {
			if (rand_double() < .2) {
				data[i] = temp;
				i++;
			}
		}
	}
}

inline double sample_quarter_circle_2_density(const vector<double> &x) {
    double area = L_PI / 4;
    double density = 1/(1-0.8*area);
    if (x[1] * x[1] + x[0] * x[0] <= 1) density /= 5;
    return density;
}


inline void sample_circle_2(int n, vector<vector<double> >& data, int d=2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < dim; j++) {
            temp[j] = rand_double();
        }
        if ((temp[1] - 0.5) * (temp[1] - 0.5) + (temp[0] - 0.5) * (temp[0] - 0.5) < 0.3 * 0.3) {
            data[i] = temp;
            i++;
        }
    }
}


inline double sample_circle_2_density(const vector<double> &x) {
    double area = L_PI * 0.3 * 0.3;
    double density = 1e-12;
    if ((x[1] - 0.5) * (x[1] - 0.5) + (x[0] - 0.5) * (x[0] - 0.5) < 0.3 * 0.3) density = 1/area;
    return density;
}

inline void sample_normal(int n, vector<vector<double> >& data, int d = 2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < dim; j++) {
            temp[j] = -1;
            while(temp[j] < 0 || temp[j] > 1) {
                temp[j] = rnorm(0.5, 0.1);
            }
        }
        data[i] = temp;
        i++;
    }
}

inline double sample_normal_density(const vector<double> &x) {
    double density = 1.0;
    for (int i = 0; i < (int)x.size(); i++) {
        density *= pnorm(x[i], 0.5, 0.1);
    }
	density /= (1-1.146606e-06);
    return density;
}

inline void sample_normal_AR1(int n, vector< vector<double> >& data, int dim, const vector<double> & normmu, double rho=0.4, double sigma=0.1, bool output=true) {
	if(normmu.size() ==dim) {
		// create covariance matrix
		Eigen::MatrixXd A(dim, dim);
		for(int i=0; i<dim; i++) {
			for(int j=0; j<dim; j++) {
				A(i,j) = sigma*sigma*exp((double)abs(i-j)* log(rho));
			}
		}
	
		// eigen decomposition of covariance matrix
		Eigen::EigenSolver<Eigen::MatrixXd> es(A);
		Eigen::MatrixXd D = es.pseudoEigenvalueMatrix();
		Eigen::MatrixXd sqrt_D = D;
		for(int i=0; i<dim; i++) {
			sqrt_D(i,i) = sqrt(D(i,i));
		}
		Eigen::MatrixXd V = es.pseudoEigenvectors();
		Eigen::MatrixXd Q(dim, dim);
		Q = V * sqrt_D * V.transpose();
		/*
		cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
		cout << "The square root eigenvalue matrix D is:" << endl << sqrt_D << endl;
		cout << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
		cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
		cout << "Finally, V * D^(1/2) * V^(-1) = " << endl << V * sqrt_D * V.transpose() << endl;
		cout << "sqrt_D * sqrt_D = "<<endl << Q * Q.transpose()<<endl;
		*/
	
		// rotate the independent data
		data.resize(n);
		int i=0;
		while(i < n) {
			data[i].resize(dim);
			// for generation of independent data
			vector<double> temp;
			temp.resize(dim);
			int flag = 0;
			while(flag < dim) {
				for(int j=0; j<dim; j++) {
					temp[j] = rnorm(0, 1);
				}
				flag=0;
				for(int j=0; j<dim; j++) {
					data[i][j] = normmu[j];
					for(int k=0; k<dim; k++) {
						data[i][j] += Q(k,j)* (temp[k]);
					}
					flag +=(int)(data[i][j] >=0 && data[i][j] <=1);
				}
			}
			i++;
		}
	
	
		if(output) {
			for(int j=0; j < dim; j++) {
				double mean=0;
				for(int i=0; i< n; i++) {
					mean += data[i][j];
				}
				mean /= n;
				cout<<"mean of dimension "<<j<<" is: "<<mean<<endl;
			}
	
			Eigen::MatrixXd B(dim, dim);
			for(int j=0; j< dim; j++) {
				for(int k=0; k< dim; k++) {
					B(j, k) = 0;
					for(int i=0; i < n; i++) {
						B(j, k) += (data[i][j]-normmu[j]) *(data[i][k]-normmu[k]);
					}	
					B(j, k) /= n;
					B(j, k) /= (sigma*sigma);
				}
			}
			cout<< "The empirical covariance matrix is:"<< endl<< B<<endl;
		}	
	}
}

inline double sample_normal_AR1_density(const vector<double> &x, const vector<double> &mu, double rho=0.4, double sigma=0.1) {
	//double rho=0.4;
	int dim = (int)x.size();
	
	// create covariance matrix
	Eigen::MatrixXd A(dim, dim);
	for(int i=0; i<dim; i++) {
		for(int j=0; j<dim; j++) {
			A(i,j) = sigma*sigma*exp((double)abs(i-j)* log(rho));
		}
	}
	double det = A.determinant();
	A = A.inverse();
	
	Eigen::VectorXd v(dim);
	for(int i=0; i < dim;i++) {
		v(i) = x[i]-mu[i];
	}
	
	Eigen::MatrixXd B(1,1);
	double density;
	density = 1/sqrt(exp(dim * log(2*M_PI))*det);
	B = v.transpose()*A*v;
	density *= exp(- B.trace() *0.5);
	
	return density;
	
}

inline void sample_mix5Dnormal(int n, vector<vector<double> >& data, int d = 5) {
    int dim = d;
    data.resize(n);
	vector<double> mu(2, 0.25);
    int i = 0;
    while (i < n) {
        vector<double> temp;
		temp.resize(dim);
		double u = rand_double();
		if(u< 0.5) {
			vector< vector<double> > temp_data;
			sample_normal_AR1(1, temp_data, 2, mu, 0.36, 0.05, false);
			temp[0] = temp_data[0][0];
			temp[1] = temp_data[0][1];
			temp[2] = -1;
			while(temp[2] <0 || temp[2]>1) {
				temp[2] = rnorm(0.25, 0.05);
			}			 
		} else{
			for(int j=0; j<3; j++) {
				temp[j]= -1;
				while(temp[j]<0 || temp[j] >1) {
					temp[j] = rnorm(0.75, 0.05);
				}
			}
		}
        for (int j = 3; j < 5 && j < dim; j++) {
            temp[j] = -1;
            while (temp[j] < 0 || temp[j] > 1) {
                temp[j] = rnorm(0.5, 0.1);
            }
        }
        for (int j = 5; j < dim; j++) {
            temp[j] = -1; 
			double temp_ind=rand_double();
			if(temp_ind < 0.5) {
				while(temp[j] <0 || temp[j]>1) {
					temp[j] = rnorm(0.35, 0.1); // for the truncated normal, partition function is 
				}			 
			}else {
				while(temp[j] <0 || temp[j]>1) {
					temp[j] = rnorm(0.6, 0.05);
				}
			}
        }
		data[i] = temp;
		i++;
	}
	return;
}

inline double sample_mix5Dnormal_density(const vector<double> &x) {
	int dim = (int)x.size();
	vector<double> mu(2, 0.25);
	vector<double> y(2, 0);
	y[0] = x[0];
	y[1] = x[1];
	double d1=1, d2=1, density=1;
	d1 = sample_normal_AR1_density(y, mu, 0.36, 0.05);
	d1 *= pnorm(x[2], 0.25, 0.05);
	for(int j=0; j<3; j++) {
		d2 *= pnorm(x[j],0.75, 0.05);
	}
	density = 0.5*d1+0.5*d2;
	for(int j=3; j<5; j++) {
		density *= pnorm(x[j], 0.5, 0.1);
	}
	for(int j=5; j<dim; j++) {
		density *= (0.5*pnorm(x[j], 0.35, 0.1) + 0.5*pnorm(x[j], 0.6, 0.05));
	}
    return density;
}



inline void sample_5normal(int n, vector<vector<double> >& data, int d=1) {
    double scale = 1;
    vector<double> sd(7,0);
    sd[0]= 0.1;
    sd[1]= 0.06;
    sd[2]= 0.044;
    sd[3]= 0.03;
    sd[4]= 0.02;
    sd[5]= 0.013;
    sd[6]= 0.01;
    double shift = 0;
    vector<double> mu(7,0);
    mu[0]= 0.5/scale+ shift;
    mu[1]= 0.4/scale+ shift;
    mu[2]= 0.33/scale+ shift;
    mu[3]= 0.28/scale+ shift;
    mu[4]= 0.26/ scale + shift;
    mu[5]= 0.24/scale+ shift;
    mu[6]= 0.21/scale+ shift;

    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp(1, -1);
        double r = rand_double();
        int choose = (int) floor(r * 7);
        while (temp[0] < 0 || temp[0] > 1) {
            temp[0] = rnorm(mu[choose], sd[choose]);
        }
        data[i] = temp;
        i++;
    }
}


inline double sample_5normal_density(const vector<double> &x) {
     double scale = 1;
    vector<double> sd(7,0);
    sd[0]= 0.1;
    sd[1]= 0.06;
    sd[2]= 0.044;
    sd[3]= 0.03;
    sd[4]= 0.02;
    sd[5]= 0.013;
    sd[6]= 0.01;
    double shift = 0;
    vector<double> mu(7,0);
    mu[0]= 0.5/scale+ shift;
    mu[1]= 0.4/scale+ shift;
    mu[2]= 0.33/scale+ shift;
    mu[3]= 0.28/scale+ shift;
    mu[4]= 0.26/ scale + shift;
    mu[5]= 0.24/scale+ shift;
    mu[6]= 0.21/scale+ shift;
    double density = 0;
    for(int k=0; k<7; k++)    density += pnorm(x[0], mu[k], sd[k])/7;
    return density;
}


inline void sample_mix2normal(int n, vector<vector<double> >& data, mixnormalparas & mixp, int d = 2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < 2; j++) {
            double u = rand_double();
            if(u<mixp.ratio){
                temp[0] = -1;
                while(temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu1, mixp.sd1);
                }

                temp[1] = -1;
                while(temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm( mixp.mu2, mixp.sd2);
                }
            }
            else{
                temp[0] = -1;
                while(temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu3, mixp.sd1);
                }
                temp[1] = -1;
                while(temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm(mixp.mu4, mixp.sd2);
                }

            }
        }
        for(int j=2; j<dim; j++)  temp[j]= rand_double();
        data[i] = temp;
        i++;
    }
}

inline double sample_mix2normal_density(const vector<double> &x, mixnormalparas & mixp) {
    double density = mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)*pnorm(x[1],mixp.mu2,mixp.sd2)+(1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1)*pnorm(x[1], mixp.mu4,mixp.sd2);
    return density;
}

//using in 'mixnormal'
inline void sample_mix2normalpro(int n, vector<vector<double> >& data, mixnormalparas & mixp, int d = 2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < 2; j++) {
            double u = rand_double();
            if (u < mixp.ratio) {
                temp[0] = -1;
                while (temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu1, mixp.sd1);
                }

                temp[1] = -1;
                while (temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm(mixp.mu2, mixp.sd2);
                }
            } else {
                temp[0] = -1;
                while (temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu3, mixp.sd1);
                }
                temp[1] = -1;
                while (temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm(mixp.mu4, mixp.sd2);
                }

            }
        }
        for (int j = 2; j < 4 && j < dim; j++) {
            temp[j] = -1;
            while (temp[j] < 0 || temp[j] > 1) {
                temp[j] = rnorm(0.5, 0.1);
            }
        }
        for (int j = 4; j < dim; j++) {
            if (rand_double() > 0.5) {
                temp[j] = -1;
                while (temp[j] < 0 || temp[j] > 1) {
                    temp[j] = rnorm(0.35, 0.1);
                }
            } else {
                temp[j] = -1;
                while (temp[j] < 0 || temp[j] > 1) {
                    temp[j] = rnorm(0.6, 0.05);
                }
            }
        }


        // the >30 dimensional case, replace the last two dimensions temp[dim-2] and temp[dim-1]
        if (dim > 31) {
            double ruf = rand_double();
            if (ruf < 0.2) {
                double a=0;
                double b=-1;
                while(a>=b){
                    a=rand_double();
                    b=rand_double();
                }
                temp[dim-2]=a;
                temp[dim-1]=b;
            }
            else if (ruf < 0.7) {
                double a=rand_double();
                double b=rand_double();
                if(a+b>=1){
                    temp[dim-2] = 0.5+0.5*a;
                    temp[dim-1] = 0.5*b;
                }
                else{
                    temp[dim-2]= 1-0.5*b;
                    temp[dim-1]= 0.5+0.5*a;
                }
            }
            else if (ruf < 0.8) {
                double a=-2;
                double b=-1;
                while(a<b){
                    a=rand_double();
                    b=rand_double();
                }
                temp[dim-2]= 1-0.5*a;
                temp[dim-1]= 0.5*b;

            }
            else {
                double a=-2;
                double b=-1;
                while(a<b){
                    a=rand_double();
                    b=rand_double();
                }
                temp[dim-2]= 0.5*a;
                temp[dim-1]= 0.5*b;
            }

        }
        data[i] = temp;
        i++;

    }
    return;

}

//using
inline double sample_mix2normalpro_density(const vector<double> &x, mixnormalparas & mixp) {
    int dim = (int)x.size();
    double density = mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)*pnorm(x[1],mixp.mu2,mixp.sd2)+(1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1)*pnorm(x[1], mixp.mu4,mixp.sd2);
    for(int j=2; j<4 && j<dim; j++)  density *= pnorm(x[j], 0.5, 0.1);
    if(dim<=31){
        for(int j=4; j<dim; j++)  density *= (0.5*pnorm(x[j],0.35,0.1)+0.5*pnorm(x[j], 0.6,0.05));
        return density;
    }
    for(int j=4; j<dim-2; j++)  density *= (0.5*pnorm(x[j],0.35,0.1)+0.5*pnorm(x[j], 0.6,0.05));
    if(x[dim-2]<x[dim-1]) density *= 0.4;
    else if((x[dim-2]+x[dim-1])>1) density *= 2;
    else if(x[dim-2]>0.5)  density *= 0.8;
    else density *= 1.6;
    return density;

}


inline double sample_mix2normalpro_sep_density(const vector<double> &x, mixnormalparas & mixp){
    double density = mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)*pnorm(x[1],mixp.mu2,mixp.sd2)+(1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1)*pnorm(x[1], mixp.mu4,mixp.sd2);
    density = density/ (mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)+ (1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1));
    density = density/ (mixp.ratio*pnorm(x[1],mixp.mu2,mixp.sd2)+ (1-mixp.ratio)*pnorm(x[1], mixp.mu4,mixp.sd2));
    return density;
}


//Confidence interval: [-5,5]            temp = (temp +5)/10;
inline void sample_skewedmixnormal(int n, vector<vector<double> >& data) {
    data.resize(n);
    for (int j = 0; j < n; j++) {
        int i = (int)floor(rand_double()*8);
        if (i == 8) i = 7;
        double temp = rnorm(3 * (pow((double)2.0 / 3.0, (int)i) - 1), pow((double)2.0 / 3.0, (int)i));
        temp = (temp +5.0)/10.0;
        data[j].push_back(temp);
        if(temp<0 && temp>1) cout<<"out of range!\n";
    }
}

inline double sample_skewedmixnormal_density(const vector<double> &x) {
    double t = x[0]*10.0-5.0;
    double density = 0;
    for(int i=0; i<8; i++)  {
        density +=   10.0/8.0 * pnorm(t, 3 * (pow(2.0 / 3.0, i) - 1),pow(2.0 / 3.0, i));         // jacobian: 10
    }

    return density;

}

inline void sample_mixstepnormal(int n, vector<vector<double> >& data, int d=6) {

    double temp;
    data.resize(n);
    for (int i = 0; i < n/2; i++) {
        data[i].push_back(rand_double() * .25 + .25);
        data[i].push_back(rand_double() * .125 + .75);
        data[i].push_back(rand_double() * .5);
        data[i].push_back(rand_double() * .5 + .5);
        temp = -1;
        while(temp < 0 || temp > 1) {
            temp = rnorm(0.3, 0.05);
        }
        data[i].push_back(temp);

        temp = -1;
        while(temp < 0 || temp > 1) {
            temp = rnorm(0.6, 0.1);
        }
        data[i].push_back(temp);
    }
    for (int i = n/2; i < n; i++) {
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
    }
}

inline double sample_mixstepnormal_density(const vector<double> &x) {
    double density1 = 0;

    if (x[0] >= 0.25 && x[0] <= 0.5 && x[1] >= 0.75 && x[1] <= 0.875 && x[2] <= 0.5 && x[3] >= 0.5) density1 += 1 / 0.25 / 0.125 / 0.5 / 0.5;
    double density2 = pnorm(x[4],0.3,0.05)*pnorm(x[5],0.6,0.1);
    return (0.5 + 0.5* density1*density2);
}

inline void sample_beta15(int n,vector<vector<double> >& data) {
    data.resize(n);
    for (int i = 0; i < n ; i++) {
        double ind = rand_double();
 //       if(ind>0.7)     data[i].push_back(rbeta(rand_double(),300,100));
 //       else  data[i].push_back(rbeta(rand_double(), 3,12));
        //if(ind>0.7)     data[i].push_back(rbeta(rand_double(),300,10));
        //else  data[i].push_back(rbeta(rand_double(), 3,12));
		int rand_seed = rand_int(0, 1073741823);
		std::mt19937_64 rng(rand_seed);
		if(ind > 0.7) data[i].push_back(rbeta_gamma_ratio(300.0, 10.0, rng));
		else data[i].push_back(rbeta_gamma_ratio(3.0, 12.0, rng));
    }
}

inline double sample_beta15_density(const vector<double> &x) {
    double density = 0.3*beta(x[0],300,10) + 0.7*beta(x[0], 3,12);
    return density;
}

inline void sample_cos(int n, vector<vector<double> >& data) {
	data.resize(n);
	int i=0;
	for(int i=0; i <n; i++) {
		bool flag = true;
		double x, y;
		while( flag) {
			x = rand_double();
			y = 2*rand_double();
			if(y <= (1+ sin(2*M_PI*x - M_PI/2))) {
				flag = false;
			}
		}
		data[i].push_back(x);
	}
}

inline double sample_cos_density(const vector<double> &x) {
	double density = 1+ sin(2*M_PI*x[0] - M_PI/2);
	return density;
}

// one-dimensional distribution for semi-circle law
inline void sample_semi_circle(int n, vector<vector<double> >& data) {
	data.resize(n);
	int i=0;
	for(int i=0; i <n; i++) {
		data[i].clear();
		bool flag = true;
		double x, y;
		while( flag) {
			x = 2* rand_double()-1;
			y = rand_double();
			if(y <= ((2.0 / M_PI) * std::sqrt( 1 - x*x)) ) {
				flag = false;
			}
		}
		data[i].push_back((x+3.0)/4.0);
	}
}

inline double sample_semi_circle_density(const vector<double> &x) {
	double density = 0.0;
	if(x[0] > 0.5 && x[0] <1) {
		density = (8.0 / M_PI) * std::sqrt( 1 - (4*x[0] - 3) * (4*x[0] -3));
	}
	return density;
}



inline void sample_1(int n,vector<vector<double> >& data) {

    data.resize(n);
    for (int i = 0; i < n / 2; i++) {
        data[i].push_back(rand_double() * .125 + 0.5);
    }
    for (int i = n / 2; i < n; i++) {
        data[i].push_back(rand_double());
    }
}

inline double sample_1_density(const vector<double> &x) {
    double density = 0.5;
    if (x[0] >= 0.5 && x[0] <= 0.625) density += 0.5 / 0.125;
    return density;
}


inline void sample_5(int n, vector<vector<double> >& data) {

    data.resize(n);
    for (int i = 0; i < n / 2; i++) {
        data[i].push_back(rand_double() * .25 + .25);
        data[i].push_back(rand_double() * .125 + .75);
        data[i].push_back(rand_double() * .5);
        data[i].push_back(rand_double() * .5 + .5);
        data[i].push_back(rand_double());
    }
    for (int i = n / 2; i < n; i++) {
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
    }
}

inline double sample_5_density(const vector<double> &x) {
    double density = 0.5;
    if (x[0] >= 0.25 && x[0] <= 0.5 && x[1] >= 0.75 && x[1] <= 0.875 && x[2] <= 0.5 && x[3] >= 0.5) density += 0.5 / 0.25 / 0.125 / 0.5 / 0.5;
    return density;
}

inline void sample_10(int n,  vector<vector<double> >& data) {
    int dim = 10;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < dim; j++) {
            temp[j] = rand_double();
        }
        if (temp[1] > 0.8 * sin(temp[0] * M_PI * 0.8)) {
            data[i] = temp;
            i++;
        } else {
            if (rand_double() < .2) {
                data[i] = temp;
                i++;
            }
        }
    }
}

inline double sample_10_density(const vector<double> &x) {
    double area = (1-cos(0.8*M_PI))/M_PI;
    double density = 1/(1-0.8*area);
    if (x[1] <= 0.8 * sin(x[0] * M_PI * 0.8)) density /= 5;
    return density;
}


// This is the function used to generate training data or test data for simulations
inline void sample_data (int n, vector< vector<double> >& data, string data_model) {
	if(data_model.find("normAR1") != string::npos) {
		size_t d_pos = data_model.find("D");
		string modelD = data_model.substr(0, d_pos);
		int data_dim = strTo<int>(modelD);
		cout<<"model dimension is: "<<modelD<<" "<<data_dim<<endl;
	
		size_t coef_pos = data_model.find(".");
		string modelCoef = data_model.substr(coef_pos);
		cout<<"autocorrelation coefficient is: "<<modelCoef<<" "<<strTo<double>(modelCoef)<<endl;
	
		vector<double> normmu(data_dim, 0.5);
	
		//unsigned seed = 283+n/100+data_dim;
		//srand(seed);
		sample_normal_AR1(n, data, data_dim, normmu, strTo<double>(modelCoef), 0.1,true);
	} else if(data_model =="mix2Dnorm") {
		mixnormalparas mixPara;
		mixPara.mu1 = 0.25;
		mixPara.mu2 = 0.25;
		mixPara.mu3 = 0.75;
		mixPara.mu4 = 0.75;
		mixPara.sd1 = 0.05;
		mixPara.sd2 = 0.05;
		mixPara.ratio = 0.4;
	
		//unsigned seed =283+n/100;
		//srand(seed);
		sample_mix2normal(n, data, mixPara, 2);
	} else if(data_model =="1Dbeta") {	
		//unsigned seed = 283+n/100;
		//srand(seed);
		sample_beta15(n, data);
	} else if (data_model =="1Dcos") {
		//unsigned seed = 283+n/100;
		//srand(seed);
		sample_cos(n, data);
	}else if (data_model == "1Dcircle") {
		sample_semi_circle(n, data);
	}else if(data_model.find("mix5Dnorm") != string::npos) {
		size_t d_pos = data_model.find_last_of("D");
		string modelD =data_model.substr(d_pos+1);
		cout<<"the total dimension is: "<<modelD<<" "<<strTo<int>(modelD)<<endl;
		int data_dim = strTo<int>(modelD);
	
		//unsigned seed = 283+n/100+data_dim;
		//srand(seed);
		sample_mix5Dnormal(n, data, data_dim);
	} else {
		cout<<"The model for data generating is not found"<<endl;
		//return 3;
	}
}


#endif //





