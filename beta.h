#ifndef BETA_H
#define BETA_H

#include "math_utils.h"

#include <random>
#include <cmath>
#include <limits>

inline double rbeta_gamma_ratio(double a, double b, std::mt19937_64& rng) {
    if (!(a > 0.0 && b > 0.0)) return std::numeric_limits<double>::quiet_NaN();

    std::gamma_distribution<double> ga(a, 1.0);
    std::gamma_distribution<double> gb(b, 1.0);

    double x = ga(rng);
    double y = gb(rng);

    // Stable computation of x/(x+y) using logs:
    // x/(x+y) = 1 / (1 + exp(log(y) - log(x))).
    double lx = std::log(x);
    double ly = std::log(y);

    double t = ly - lx;
    if (t > 0) {
        // avoid overflow of exp(t) when t is large
        double e = std::exp(-t);
        return e / (1.0 + e);
    } else {
        double e = std::exp(t);
        return 1.0 / (1.0 + e);
    }
}

static const double IntDis = 0.00001;


inline double  FBeta(double a, double b){
	return exp(lgamma_c(a)+lgamma_c(b)-lgamma_c(a+b));
}

//beta(a,b) function of x
inline double beta(double x, double a, double b){
    if(x<0 || x>1) return 0;
    return(pow(x,a-1)*pow((1-x), b-1))/FBeta(a,b);
}

inline double cdfbeta(double left, double right, double a, double b){
	double cdf =0;
	for(double x=left;x<right; x+=IntDis){
		cdf = cdf + beta(x, a,b)*IntDis;
	}
	return(cdf);
}


inline double rbeta(double u, double a, double b){
	double r = 0;
	double cdf=0;
	while(cdf < u && r<=1){
            r = r+IntDis;
            cdf = cdf + beta(r-0.5*IntDis,a,b)*IntDis;

	}
        if(cdf<u)  cerr<<"wbs!: u="<<u<<" cdf="<<cdf<<"\n";
	return (r);

}

#endif

