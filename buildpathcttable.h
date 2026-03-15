/*
 * File:   buildpathcttable.h
 * Author: Mouse
 *
 * Created on October 27, 2011, 3:52 PM
 */

#ifndef BUILDPATHCTTABLE_H
#define	BUILDPATHCTTABLE_H

#include "sampling.h"
#include "output.h"



// given the dimension of the first cut and the patterns for both sides of the cut, add a partition in the map together with its score
//struct shrink_reg{
//	vector<usint_mask> reg_code;
//	int num;
//};
inline void  add(map< OnePartition , double, CompairSRegs>& mpathctmap, int d, OnePartition P1, double lnumpath1, OnePartition P2, double lnumpath2 ){
    OnePartition P;
    for(int i=0; i<(int)P1.sregs.size(); i++){
        shrink_reg sr =  P1.sregs[i];
        sr.reg_code[d].mask  = (sr.reg_code[d].mask <<1) +1;
        P.sregs.push_back(sr);

    }
    for(int j=0; j<(int)P2.sregs.size(); j++){
        shrink_reg sr =  P2.sregs[j];
        sr.reg_code[d].x = sr.reg_code[d].mask +1+sr.reg_code[d].x;
        sr.reg_code[d].mask  = (sr.reg_code[d].mask <<1) +1;
        P.sregs.push_back(sr);

    }
    double lnumpath = lnumpath1 + lnumpath2 + stirling((int)P1.sregs.size()-1 + (int)P2.sregs.size()-1)-stirling((int)P1.sregs.size()-1)-stirling((int)P2.sregs.size()-1) ;

    sort(P.sregs.begin(), P.sregs.end(), myfunction3);

    map<OnePartition , double, CompairSRegs>::iterator it;
    it= mpathctmap.find(P);
    if(it==mpathctmap.end()){
        mpathctmap.insert(pair<OnePartition, double>( P, lnumpath));

    }
    else {
        it->second = log(exp(it->second) + exp(lnumpath));

    }

}


// given the first cut dimension and the number of cuts in each regions.
inline void  add(map< OnePartition , double, CompairSRegs>& mpathctmap, int d, map< OnePartition , double, CompairSRegs> map1, map< OnePartition , double, CompairSRegs> map2){
    map<OnePartition , double, CompairSRegs>::iterator it1, it2;
    for ( it1 = map1.begin() ; it1 != map1.end(); it1++ ){
        for ( it2 = map2.begin() ; it2 != map2.end(); it2++ ){
            add(mpathctmap, d, it1->first, it1->second, it2->first,it2->second );
        }
    }
}

// m is the number of cuts
inline void onebuildpathcttable(int m, int dim, vector< map< OnePartition , double, CompairSRegs> >& pathctmap ){
    if(m==0){            // add the samplespace together with the score
        OnePartition P;
        shrink_reg sr;
        usint_mask um;
        um.x= 0;
        um.mask = 0;
        for(int d=0; d<dim; d++)   {
            sr.reg_code.push_back(um);
        }
        P.sregs.push_back(sr);
        map< OnePartition , double, CompairSRegs> map0;
        map0.insert(pair<OnePartition , double> (P, 0.0));
        pathctmap.push_back(map0);
        return;

    }

    int firstpartcuts, secondpartcuts;
    map< OnePartition , double, CompairSRegs> mpathctmap;
    mpathctmap.clear();
    for(firstpartcuts=0; firstpartcuts< m; firstpartcuts++ ){
        secondpartcuts = m -1- firstpartcuts;
        for(int d=0; d<dim; d++){
            add(mpathctmap, d, pathctmap[firstpartcuts], pathctmap[secondpartcuts]);
        }

    }
    pathctmap.push_back(mpathctmap);
    return;

}

inline void buildpathcttable(int m, int dim, vector< map< OnePartition, double, CompairSRegs> >& pathctmap ){
    pathctmap.clear();
    for(int i=0; i<m; i++){
         cerr<<"add="<<i<<endl;
    //    cout<<"i="<<i<<"\n";
        onebuildpathcttable(i,dim, pathctmap);
       cout<<"pathctmap.size()="<<(int)pathctmap.size()<<"\n";

    }

}



inline void printout( map< OnePartition , double, CompairSRegs >& pathctmap){
    map<OnePartition , double, CompairSRegs>::iterator it;
    int s=0;
    for( it= pathctmap.begin(); it!=pathctmap.end(); it++){
        print_partition(it->first);
        s++;
       cout<<"lpathct="<<it->second<<"\n\n";
    }
    cout<<"totalnum="<<s<<"\n";

}




#endif	/* BUILDPATHCTTABLE_H */

