#include "tree.h"
#include "sampling.h"
#include "timer_ed520.h"
#include <iomanip>
using namespace std;
extern TIMER timer;

bool compare_weight ( int i , int j ){
    return ( i<j );
}
bool operator==(  const vector<int> &a , const vector<int> &b) {
    if( a.size() != b.size()){
        return false;
    }
    for(unsigned i=0 ; i < a.size() ; i++){
        if(a[i] != b[i])
            return false;
    }
    return true;
}

Partition::Partition(){
    w[0]=0;
    w[1]=0;
    score[0]=0.0;
    score[1]=0.0;
}

double tree::lBPscore_t(Partition & P){
    //Update BPscore for partition P by adding scores of two new region and minusing the score of the cut region.
    int cutreg = P.cutregs[List_index];
    P.score[List_index] += lgamma( 0.5 + P.region[List_index][cutreg]->length);
    P.score[List_index] += lgamma( 0.5 + P.region[List_index][ P.region[List_index].size()-1]->length);
    P.score[List_index] -= lgamma( 0.5 + (P.region_replace->length));
    
    P.score[List_index] -= ((P.region[List_index][cutreg])->length)*lprod_usint_mask(P.region[List_index][cutreg]->Region_mask );
    P.score[List_index] -= ((P.region[List_index][P.region[List_index].size()-1])->length)*lprod_usint_mask(P.region[List_index][P.region[List_index].size()-1]->Region_mask );
    P.score[List_index] += (P.region_replace->length)*lprod_usint_mask( P.region_replace->Region_mask);
    
    return  P.score[List_index];
}

double tree::getmaxlBP_t(  vector<Partition> &Pts){

    double maxlBP = lBPscore_t(Pts[0]);
    modelBP=0;
    double temp;
    for (unsigned i = 1; i< Pts.size(); i++) {
        temp= lBPscore_t(Pts[i]);
        if(temp>maxlBP){
            maxlBP = temp;
            modelBP=i;
        }
    }
    return (maxlBP);

}

void tree::tree_nonpara_init(){
#ifdef PROFILE
    timer.TimerStart_cpu(timer.sTime_init_NC);
#endif

    beta =20;
    maxpercentage =0.9;
    Samplesize = Data.size(); 
    // resample parameter 
    Steps=4; 
    Resampling=2;
    // end  
    hit_in_count=0;
    Data_array= new double[Dim*Samplesize];
    Point_index= new int[Samplesize*n*SCALE];
    
    Last_point_index=Samplesize;
    Tmp_count= new int [Dim];
    Tmp_unique = new bool [Dim];
    //Root.sub_length.clear();
    // Root initialzation
    //Bsp_Tree.clear();
    //Root.Node_child.clear();
    Root.um.mask=0;
    Root.um.x=0;
    for(int i=0;i<Dim;i++)
    {
        Root.Region_mask.push_back(Root.um);  // looks fine for the region boundary
    }
    Root.Data_point_start = 0;
    Root.c_original.resize(Dim,0);
    Root.c_modified.resize(Dim,0);
    Root.unique_dim = new bool [Dim];
    Root.length = Samplesize;
    Root.Node_parent = NULL;
    //end root init
    //Bsp_Tree.clear();
    Bsp_Tree.push_back(&Root);
    for(int i=0;i<Samplesize;i++){
        for(int j=0;j<Dim;j++){
            Data_array[i*Dim+j]= Data[i][j];
        }
    }
    for(int i=0;i<Samplesize;i++){
        Point_index[i]= Dim*i;
    }
    //Pts.clear();
    double init_score=lgamma( 0.5 +Root.length );
    List_index=0;
    List_index_best=0;
    Partition partits[n]; //declare n partitions
    for(int i=0;i<n;i++){
        partits[i].first=0;
        partits[i].region_1.push_back(&Root);
        partits[i].region_2.push_back(NULL);
        partits[i].region.push_back(partits[i].region_1);
        partits[i].region.push_back(partits[i].region_2);
        partits[i].PartitionsID=i;
        partits[i].score[0] = init_score;
        Pts.push_back(partits[i]);
    }
    for(int i=0 ; i<n ; i++ )
        Root.Partitions.push_back(&Pts[i]);
    
    subtree_root.push_back(&Root);

#ifdef PROFILE
    timer.TimerFinish_cpu(timer.tTime_init_NC,timer.sTime_init_NC);
#endif
}
void tree::para_share_init(){

    Dim=1;
    Point_index= new int[n*Samplesize];
    
    Data_array = new double[Samplesize];

    Tmp_count= new int [Dim];
    Tmp_unique = new bool [Dim];

}

void tree::tree_para_init( ){
#ifdef PROFILE
    timer.TimerStart_cpu(timer.sTime_init_C);
#endif

    Root.sub_length.clear();
    Last_point_index=Samplesize;
    Root.Node_child.clear();    
    Bsp_Tree.clear();
    //  initializate boundary for each dim

    // flag area
    hit_in_count=0;

    Root.Region_mask.clear();
    Root.um.mask=0;
    Root.um.x=0;
    for(int i=0 ; i < Dim ; ++i){
        Root.Region_mask.push_back(Root.um);
    }
    Root.unique_dim = new bool[Dim];
    Root.Data_point_start=0;
    Root.length = Samplesize;
    Root.c_original.resize(Dim,0);
    Root.c_modified.resize(Dim,0);
    Bsp_Tree.clear();
    Bsp_Tree.push_back(&Root);
    double init_score=lgamma(0.5 + Root.length);
    List_index=0;
    List_index_best=0;
    Pts.clear();
    Partition para_parts[n];
    for(int i=0;i<n;i++){
        para_parts[i].first=0;
        para_parts[i].region_1.clear();
        para_parts[i].region_1.push_back(&Root);
        para_parts[i].region_2.clear();
        para_parts[i].region_2.push_back(NULL);
        para_parts[i].region.clear();
        para_parts[i].region.push_back(para_parts[i].region_1);
        para_parts[i].region.push_back(para_parts[i].region_2);
        para_parts[i].PartitionsID=i;
        para_parts[i].score[0] = init_score;
        Pts.push_back(para_parts[i]);
    }

#ifdef PROFILE
    timer.TimerFinish_cpu(timer.tTime_init_C,timer.sTime_init_C);
#endif
}

void Region_Node::destroy(){
    Partitions.clear();
    delete [] unique_dim;
    c_original.clear();
    c_modified.clear();
    c_local_max.clear();
    sub_length.clear();
    Dim_has_choose.clear();
    Node_child.clear();
    Region_mask.clear();
}

int Region_Node::cut_dim(){   		//+ 412
    int cut_dimension = -1;
    if ( Node_child.size() != 0){
        for (unsigned int d=0; d<Region_mask.size(); d++){
            pair<double,double> subreg0_range = convert_ranges(Node_child[0]->Region_mask[d]);    
            pair<double,double> subreg1_range = convert_ranges(Node_child[1]->Region_mask[d]);
            if (subreg0_range.first != subreg1_range.first){
                cut_dimension = d;
                break;
            }
        }
    }
    return cut_dimension;
}

double Region_Node::find_volume(){	//+ 412
    double vol = 1;
    for (unsigned int d=0; d<Region_mask.size(); d++){
        pair<double,double> range = convert_ranges(Region_mask[d]);
        vol *= (range.second - range.first);
    }
    return vol;
}

void tree::tree_para_destroy(){
    delete [] Point_index;
    delete [] Data_array;
    delete [] Tmp_count;
    delete [] Tmp_unique;
}

void tree::SIS_tree( ){
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerStart_cpu(timer.sTime_C_SIS);
    else
        timer.TimerStart_cpu(timer.sTime_NC_SIS);
#endif

    beta =1;
    cout<<"beta = "<<beta<<endl;
    modelBP=0;
    currlBPscore=-1.0;
    P_bestlevel = Pts[modelBP];     //added by thchiu, need the case when you don't cut 
    maxlBPscore=currlBPscore;
    maxlevel=0;
    cout<<"Dim "<<Dim<<endl;
    for(int i=0;i<Levels;i++){
        cout<<"\n\n#i=" <<i<<endl;
        if((i>0)&&(Resampling)&&(i%Steps==0)  /* 's times and resample mechanism is enhabled   */  ){   
            Resample_function();
        }
        bool exhaust = Sample_tree();
        if(exhaust){
            break;  
        }
#ifdef PROFILE
        timer.TimerStart_cpu(timer.sTime_BPscore);
#endif
        currlBPscore = -lgamma_c(Samplesize + (i + 2) / 2.0)- (i + 2) * lgamma_c(0.5) + lgamma_c(0.5 * (i + 2)) - beta * (i + 2);

        currlBPscore = currlBPscore + getmaxlBP_t(Pts);

        cout<<currlBPscore << "\n";
        
        if (currlBPscore > maxlBPscore) {
            maxlBPscore = currlBPscore;
            P_bestlevel = Pts[modelBP]; // seems OK! 
            List_index_best = List_index; 
            cout<<"modelBP="<<modelBP<<endl;
            maxlevel=i;
            cout<<"newmax level="<<i+2<<endl;
        }
#ifdef PROFILE
        timer.TimerFinish_cpu(timer.tTime_BPscore,timer.sTime_BPscore);
#endif
        //Add is_milestone, thchiu
/*
        if( is_milestone(i, p.levels) && temp_ofilename!="" ) {
            string ofilename2 =temp_ofilename+"_level_"+toStr<int>(i)+".txt";
            ofstream outfile (ofilename2.c_str());
            print_partition(outfile, P_bestlevel, temp_mmax, temp_mmin);
        }
*/
        if(i>maxlevel+10) {
            cout<<"no improvement in 10 levels!\n";
            break;
        }
    }

#ifdef PROFILE
    if(Dim == 1)
        timer.TimerFinish_cpu(timer.tTime_C_SIS, timer.sTime_C_SIS);
    else
        timer.TimerFinish_cpu(timer.tTime_NC_SIS, timer.sTime_NC_SIS);
#endif
}


bool tree::is_milestone(int i , int plevels ){

    if(i==(int)((double)plevels/10.0)) return true;
    if(i==(int)((double)plevels/5.0)) return true;
    if(i==(int)((double)plevels/3.0)) return true;
    if(i==(int)((double)plevels/2.0)) return true;
    if(i==(int)((double)plevels*2.0/3.0)) return true;
    if(i==(int)((double)plevels*4.0/5.0)) return true;
    return false;

}

bool tree::Sample_tree(){
    smooth =0.0; 
    for(unsigned nIter = 0  ; nIter < Pts.size(); ++nIter ){
        
        //cout<<"partition "<<nIter<<endl;
        get_info_from_onetree(nIter);    

        if( Pts[nIter].cmax < 1- MAXNUM ){
            return true;
        }
        rand_choose_one_nextleveltree( nIter );    
    }

    double largestwt = Pts[0].w[List_index];
    for( unsigned nIter=1; nIter< Pts.size(); nIter++ ){
        if(Pts[nIter].w[List_index]>largestwt)  largestwt= Pts[nIter].w[List_index];
    }
    for (unsigned nIter = 0; nIter< Pts.size(); nIter++) {  
        Pts[nIter].w[List_index] = Pts[nIter].w[List_index] - largestwt + 5;
    }

    return false;   
}


void tree::get_info_from_onetree(int nIter){
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerStart_cpu(timer.sTime_get_info_C);
    else
        timer.TimerStart_cpu(timer.sTime_get_info_NC);
#endif
    if( Pts[nIter].first ==0 ){
        Count(Pts[nIter].region[List_index][0]);  
        Calculate_weight(nIter,0);                   
        Pts[nIter].first =1;
    }
    else {
        unsigned Last_one_reg = Pts[nIter].region[List_index].size()-1;
        //  First one Region
        Count(Pts[nIter].region[List_index][Pts[nIter].cutregs[List_index]]);
        Calculate_weight (nIter,Pts[nIter].cutregs[List_index]); 
        //  Second one Region 
        Count(Pts[nIter].region[List_index][Last_one_reg]);
        Calculate_weight(nIter,Last_one_reg);
    }
    Two_Maxc_in_one_part(Pts[nIter]);

#ifdef PROFILE
    if(Dim == 1)
        timer.TimerFinish_cpu(timer.tTime_get_info_C , timer.sTime_get_info_C);
    else
        timer.TimerFinish_cpu(timer.tTime_get_info_NC , timer.sTime_get_info_NC);
#endif
}


//now we implement the discrete is false first we will lately to implement a true version !!! 
void tree::Count(   Region_Node * region_cut  ){
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerStart_cpu(timer.sTime_Count_C);
    else
        timer.TimerStart_cpu(timer.sTime_Count_NC_cpu);
#endif

    //case 1 : counting has done before
    if(region_cut->sub_length.size()!=0){
#ifdef PROFILE
        if(Dim != 1) timer.HIT_COUNT++;
#endif
        for(int iter =0 ; iter < Dim ; iter ++  ){
            Tmp_count[iter] = region_cut->sub_length[iter];   // direct access
            Tmp_unique[iter]= region_cut->unique_dim[iter]; 
        }
        hit_in_count=1;
    }
    //case 2 : count here !
    else {   
#ifdef PROFILE
        if(Dim != 1) timer.CPU_COUNT++;
#endif
        double * Tmp_upbound;
        Tmp_upbound = new double [Dim];
        usint_mask Tmp_reg_code;
        for(int d=0;d<Dim;d++){
            Tmp_count[d] =0;
            Tmp_reg_code.x=(region_cut->Region_mask[d].x<<1 );
            Tmp_reg_code.mask=(region_cut->Region_mask[d].mask<<1) +1;
            pair<double ,double> ranges =convert_ranges(Tmp_reg_code);
            Tmp_upbound[d]=ranges.second;
            if(discrete) Tmp_unique[d] =true;
            else         Tmp_unique[d] =false;
        }
        for (int i=0;i<region_cut->length;i++){
            for(int j=0;j<Dim;j++){
                if(Data_array[ Point_index[region_cut->Data_point_start+i]+j] <= Tmp_upbound[j]  ){
                    ++Tmp_count[j];
                }

                if(Tmp_unique[j] && fabs(Data_array[ Point_index[region_cut->Data_point_start]+j]-Data_array[ Point_index[region_cut->Data_point_start+i]+j]) > 1e-10){
                    Tmp_unique[j] = (Tmp_unique[j] && false);
                }
            }
        }
        for(int i=0;i<Dim;i++){
            region_cut->sub_length.push_back(Tmp_count[i]); 
            region_cut->unique_dim[i]=Tmp_unique[i];
        }
        hit_in_count=0;
    }

#ifdef PROFILE
    if(Dim ==1)
        timer.TimerFinish_cpu(timer.tTime_Count_C,timer.sTime_Count_C);
    else
        timer.TimerFinish_cpu(timer.tTime_Count_NC_cpu,timer.sTime_Count_NC_cpu);
#endif
}

void tree::Region_Create(Region_Node *region_cut,Region_Node *&New_Region_L,Region_Node *&New_Region_R,int i){

    ///Left child
    New_Region_L = new Region_Node ;
    New_Region_L->c_original.resize(Dim,0);  // necessary to be implemented in the count function  
    New_Region_L->c_modified.resize(Dim,0);  // same as the up !!!
    New_Region_L->unique_dim= new bool [Dim];
    //assign boundary
    New_Region_L->Region_mask=region_cut->Region_mask;  // assign first    
    New_Region_L->Region_mask[i].x = New_Region_L->Region_mask[i].x<<1;
    New_Region_L->Region_mask[i].mask=(New_Region_L->Region_mask[i].mask<<1)+1;
    New_Region_L->length= region_cut->sub_length[i];    // Tmp_count[i];   // assign data amount

    ///Right child
    New_Region_R = new Region_Node ;
    New_Region_R->c_original.resize(Dim,0);  // important! fine it can be in this function !! 
    New_Region_R->c_modified.resize(Dim,0);  // important! as the above  
    New_Region_R->unique_dim= new bool [Dim];
    //assign boundary
    New_Region_R->Region_mask=region_cut->Region_mask;
    New_Region_R->Region_mask[i].x = New_Region_R->Region_mask[i].x<<1;
    New_Region_R->Region_mask[i].mask=(New_Region_R->Region_mask[i].mask<<1) +1;
    New_Region_R->Region_mask[i].x++;              
    New_Region_R->length=region_cut->length-region_cut->sub_length[i];     // assign data amount

    //update region
    region_cut->Node_child.push_back(New_Region_L);
    region_cut->Node_child.push_back(New_Region_R);
    New_Region_L->Node_parent = region_cut;
    New_Region_R->Node_parent = region_cut;
    //update the Bsp_Tree
    Bsp_Tree.push_back( New_Region_L);
    Bsp_Tree.push_back( New_Region_R);
}

void tree::Two_Maxc_in_one_part( Partition & P  ){
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerStart_cpu(timer.sTime_Two_Maxc_in_one_part_C);
    else
        timer.TimerStart_cpu(timer.sTime_Two_Maxc_in_one_part_NC);
#endif

    unsigned region_size = P.region[List_index].size();
    
    // new method here, record the maximum weight & 2nd maximum weight
    P.cmax  = -2*MAXNUM;
    P.cmax2 = -2*MAXNUM;
    int max_index=-1;

    for(unsigned i=0; i < region_size ; i++){
        if(P.region[List_index][i]->c_local_max[0] > P.cmax){
            P.cmax2 = P.cmax;
            P.cmax = P.region[List_index][i]->c_local_max[0];
            max_index = i;
        }
        else if(P.region[List_index][i]->c_local_max[0] < P.cmax && P.region[List_index][i]->c_local_max[0] > P.cmax2){
            P.cmax2= P.region[List_index][i]->c_local_max[0];
        }
    }
    if(Dim != 1){
        if( P.cmax2 < P.region[List_index][max_index]->c_local_max[1]){
            P.cmax2 = P.region[List_index][max_index]->c_local_max[1];
        }
    }
    //new method finished
    
    if(P.cmax < 1 - MAXNUM)
    {
        cout<<"no regions no dimensions to split !!" <<endl;
        return; 
    }

    if( P.cmax2 < 1-MAXNUM ) {
        return;
    }

#ifdef DEBUG
    for(unsigned i = 0 ; i < region_size ; i ++ ){
        for(int j=0 ; j <Dim ; j++) {
            if(isnan(  P.region[List_index][i]->c_original[j] )){
                cout<<"region="<<i<<"dim="<<j<<"c_original["<<j<<"]=NAN at a"<<endl;
            }
        }
    }
#endif

    if( maxpercentage < 0  ){
        smooth=20;
    }
    else{
        smooth = (P.cmax-P.cmax2)/ log(maxpercentage / (1 - maxpercentage));
    }
    if(isinf(smooth)){ 
        cout<<"smooth==inf! cmax="<<P.cmax<<" cmax2="<<P.cmax2<<endl;
    }


    if (smooth > 1) {
        for (unsigned i = 0; i < region_size; i++) {
            for(int j = 0 ; j <Dim ; j++  ){
                if( P.region[List_index][i]->c_original[j] > 1-MAXNUM )  {
                    P.region[List_index][i]->c_modified[j] = P.region[List_index][i]->c_original[j]/smooth;
                }
                else
                    P.region[List_index][i]->c_modified[j] = P.region[List_index][i]->c_original[j];
            }
        }
    }
    else{
        for (unsigned i = 0; i < region_size; i++) {
            for(int j = 0 ; j <Dim ; j++  ){
                P.region[List_index][i]->c_modified[j] = P.region[List_index][i]->c_original[j] ;
            }
        }
    }

#ifdef DEBUG 
    // checking whether it is NAN or not 
    for(unsigned i = 0; i < region_size; i++) {
        for(int j=0 ; j< Dim ; j++){
            if(isnan( P.region[List_index][i]->c_modified[j] )) {
                cout<<"region="<<i<<"dim="<<j<<"c_modified["<<j<<"]=NAN at a"<<endl; 
                cout<<"smooth="<<smooth<<endl;
            }
        }
    }
    // finish the second checking !!! 
#endif

#ifdef PROFILE
    if(Dim == 1)
        timer.TimerFinish_cpu(timer.tTime_Two_Maxc_in_one_part_C,timer.sTime_Two_Maxc_in_one_part_C);
    else
        timer.TimerFinish_cpu(timer.tTime_Two_Maxc_in_one_part_NC,timer.sTime_Two_Maxc_in_one_part_NC);
#endif
}


void tree::Calculate_weight( int  nIter ,  int reg_index  ){
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerStart_cpu(timer.sTime_Calculate_weight_C);
    else
        timer.TimerStart_cpu(timer.sTime_Calculate_weight_NC);
#endif

    if( hit_in_count ==0 ){
        int samplesize = Pts[nIter].region[List_index][reg_index]->length;//the number in this mother region!
        for(int i= 0 ; i <Dim ; ++i ){
            Pts[nIter].region[List_index][reg_index]->c_original[i]=0;
            double beta2 =0.5;
            if( (double) samplesize / 200.0 < 0.5  ){
                beta2 = max( 0.1 , (double) samplesize /200.0  ); // heuristuc pseudo ct
            }
            if(Tmp_unique[i]){
                Pts[nIter].region[List_index][reg_index]->c_original[i] = -MAXNUM;
            }
            else{
                Pts[nIter].region[List_index][reg_index]->c_original[i] += lgamma_c(beta2 + Tmp_count[i]  );
                Pts[nIter].region[List_index][reg_index]->c_original[i] += lgamma_c(beta2 + (samplesize - Tmp_count[i]) );

            }
            if( Pts[nIter].region[List_index][reg_index]->c_original[i]  > 1 - MAXNUM  ){
                Pts[nIter].region[List_index][reg_index]->c_original[i] -=lgamma_c(beta2 + samplesize);
                Pts[nIter].region[List_index][reg_index]->c_original[i] +=samplesize*log(2.0);
            }
        }

        ///New method, record the maximum c_original & 2nd maximum c_original
        if(Dim == 1){  //There are only 1 dimension, so only one element in c_original
            Pts[nIter].region[List_index][reg_index]->c_local_max.resize(1,0);
            Pts[nIter].region[List_index][reg_index]->c_local_max[0] = Pts[nIter].region[List_index][reg_index]->c_original[0];
        }
        else{
            double temp_max= -MAXNUM;
            double temp_max2= -2*MAXNUM;

            // looking for the Max and second max
            for(unsigned  i=0; i<Pts[nIter].region[List_index][reg_index]->c_original.size(); i++){
                if(Pts[nIter].region[List_index][reg_index]->c_original[i] > temp_max ){
                    temp_max2 = temp_max;
                    temp_max = Pts[nIter].region[List_index][reg_index]->c_original[i];
                }
                else if(Pts[nIter].region[List_index][reg_index]->c_original[i] < temp_max  &&  Pts[nIter].region[List_index][reg_index]->c_original[i] > temp_max2 ){
                    temp_max2 = Pts[nIter].region[List_index][reg_index]->c_original[i] ;
                }
            }
            Pts[nIter].region[List_index][reg_index]->c_local_max.resize(2,0);
            Pts[nIter].region[List_index][reg_index]->c_local_max[0] = temp_max;
            Pts[nIter].region[List_index][reg_index]->c_local_max[1] = temp_max2;
        }
    }

#ifdef PROFILE
    if(Dim == 1)
        timer.TimerFinish_cpu(timer.tTime_Calculate_weight_C,timer.sTime_Calculate_weight_C);
    else
        timer.TimerFinish_cpu(timer.tTime_Calculate_weight_NC,timer.sTime_Calculate_weight_NC);
#endif
}

void tree::rand_choose_one_nextleveltree(int nIter){
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerStart_cpu(timer.sTime_rand_choose_C);
    else
        timer.TimerStart_cpu(timer.sTime_rand_choose_NC);
#endif

    vector <double> ec( Dim*Pts[nIter].region[List_index].size(),0);
    unsigned region_size= Pts[nIter].region[List_index].size();
    //find c_max of all regions
    double rand_cmax= Pts[nIter].region[List_index][0]->c_modified[0];
    for(unsigned i=0;i<region_size;i++){ //must begin at region 0, because we need all the c_modified within all the dimension in region 0 
        for(int j=0;j<Dim;j++){
            if(Pts[nIter].region[List_index][i]->c_modified[j]>rand_cmax){
                rand_cmax = Pts[nIter].region[List_index][i]->c_modified[j];
            }
        }
    }
    double smallerthanm=0;
    for(unsigned i=0;i<region_size;i++){
        for(int j=0;j<Dim;j++){
            if( Pts[nIter].region[List_index][i]->c_modified[j]<1-MAXNUM){
                ec[i*Dim+j]=0;
            }
            else{
                //ec[i*Dim+j]= exp((Pts[nIter].region[List_index][i]->c_modified[j])-rand_cmax+5);
                ec[i*Dim+j] = exp((Pts[nIter].region[List_index][i]->c_modified[j])-rand_cmax);
                smallerthanm += ec[i*Dim+j];
            }
            if(isnan(Pts[nIter].region[List_index][i]->c_modified[j])) 
                cout<< " is.nan(ec[" << i << "]), Pts[nIter].region[List_index][i]->c_original[j]=" <<Pts[nIter].region[List_index][i]->c_modified[j]<<endl;
        }
    }
    int &cutreg = Pts[nIter].cutregs[List_index];
    int cut = rand_int(ec);       // randomly choose a cut 
    int cutvar = cut%Dim ;        // which dim to cut 
    cutreg= (cut-cutvar)/Dim;     // which region to cut! the new region is alive       

    //need to modulation
    double c_cut=0.0 ;
    c_cut=Pts[nIter].region[List_index][cutreg]->c_modified[cutvar];
    double logsum_c ; 
    
    logsum_c = rand_cmax + log(smallerthanm); 
    if(smooth>1){
        Pts[nIter].w[List_index] = Pts[nIter].w[List_index] + (smooth -1) *c_cut + logsum_c;
    }
    else{
        Pts[nIter].w[List_index]= Pts[nIter].w[List_index] + logsum_c;
    }
    //end modulation    
    Region_Node * region_select = Pts[nIter].region[List_index][cutreg];
    Pts[nIter].region_replace = region_select; //store the replaced region
    
    double Tmp_upbound;
    usint_mask mytemp;
    mytemp.x=(region_select->Region_mask[cutvar].x<<1);
    mytemp.mask =(region_select->Region_mask[cutvar].mask<<1)+1;
    pair <double,double> ranges = convert_ranges(mytemp);
    Tmp_upbound = ranges.second;

    //  new method begins here 
    if(region_select-> Dim_has_choose.size() == 0 ){  // here means the region we select is raw !!!
        int data_head =  region_select->Data_point_start;   
        int data_tail =  data_head+region_select->length;   
        int Tmp_count=0;
        for(;data_head != data_tail;){
            if(Data_array[ Point_index[data_head]+cutvar] > Tmp_upbound){
                --data_tail;
                int tmp= Point_index[data_head];
                Point_index[data_head]=Point_index[data_tail];
                Point_index[data_tail]=tmp;

            }
            else{
                Tmp_count++;
                ++data_head;
            }
        }

        if(data_head <region_select->Data_point_start + region_select->length)
            if(Data_array[ Point_index[data_head]+cutvar] <= Tmp_upbound)
                Tmp_count++;
        
        Region_Node *New_Region_L=NULL;
        Region_Node *New_Region_R=NULL;
        Region_Create( region_select , New_Region_L , New_Region_R, cutvar);
        unsigned Last_two = region_select->Node_child.size()-2;
        unsigned Last_one = region_select->Node_child.size()-1;  
        region_select->Node_child[Last_two]->Data_point_start = region_select->Data_point_start;
        region_select->Node_child[Last_one]->Data_point_start = region_select->Data_point_start+Tmp_count;
        Pts[nIter].region[List_index][cutreg]=region_select->Node_child[Last_two];
        Pts[nIter].region[List_index].push_back(region_select->Node_child[Last_one]);
        region_select->Dim_has_choose.push_back(cutvar);
    }
    else {
        int test_t=0; 
        for(unsigned i=0;i< region_select->Dim_has_choose.size();i++){
            if(cutvar==region_select->Dim_has_choose[i]){
                test_t=1;  // congrdulation !! the region exist already !!
                Pts[nIter].region[List_index][cutreg] = region_select->Node_child[2*i];
                Pts[nIter].region[List_index].push_back(region_select->Node_child[2*i+1]);
                break;
            }
        }
        if (test_t ==0 ){  // the worst case happen 
            if(Last_point_index+region_select->length>n*Samplesize*SCALE)  // the max size we can afforad !
            {
                cout<<"Memory space is not enough!!"<<endl;
                cout<<"Last_point_index is "<<Last_point_index<<endl;
                Data_point_defragment();
                cout<<"defragment Last_point_index is "<<Last_point_index<<endl;
            }
            for(int i=0;i<region_select->length;i++){   //region index  copy
                Point_index[ Last_point_index + i]=Point_index[region_select->Data_point_start+i];// the index is copied now
            }

            int data_head = Last_point_index;
            int data_tail = Last_point_index+region_select->length;
            int Tmp_count =0;
            for( ; data_head!= data_tail ;  ){
                if(Data_array[Point_index[data_head]+ cutvar]>Tmp_upbound){
                    --data_tail;
                    int tmp = Point_index[data_head];
                    Point_index[data_head]= Point_index[data_tail];
                    Point_index[data_tail]=tmp;
                }
                else{
                    Tmp_count++;
                    ++data_head;
                }
            }
            if(data_head <Last_point_index + region_select->length)
                if( Data_array[Point_index[data_head]]+ cutvar<= Tmp_upbound)
                    Tmp_count++;
            
            Region_Node *New_Region_L=NULL;
            Region_Node *New_Region_R=NULL;
            Region_Create( region_select , New_Region_L , New_Region_R, cutvar);
            unsigned Last_two = region_select->Node_child.size()-2;
            unsigned Last_one = region_select->Node_child.size()-1;
            region_select->Node_child[Last_two]->Data_point_start = Last_point_index;
            region_select->Node_child[Last_one]->Data_point_start = Last_point_index+Tmp_count;
            Pts[nIter].region[List_index][cutreg]= region_select->Node_child[Last_two] ;
            Pts[nIter].region[List_index].push_back(region_select->Node_child[Last_one]);
            Last_point_index += region_select->length;
            region_select->Dim_has_choose.push_back(cutvar);
            
            if(Dim != 1){
                subtree_root.push_back(region_select->Node_child[Last_two]);
                subtree_root.push_back(region_select->Node_child[Last_one]);
            }
        }
    }
    
    if(Dim != 1){
        Erase_partition( region_select , nIter );
        Pts[nIter].region[List_index][cutreg]->Partitions.push_back( &Pts[nIter]  );
        int last_region = Pts[nIter].region[List_index].size()-1;
        Pts[nIter].region[List_index][last_region]->Partitions.push_back(&Pts[nIter]);
    }
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerFinish_cpu(timer.tTime_rand_choose_C , timer.sTime_rand_choose_C);
    else
        timer.TimerFinish_cpu(timer.tTime_rand_choose_NC ,timer.sTime_rand_choose_NC);
#endif
}

void tree::Data_point_defragment(){
    stack <Region_Node *>  root_stack; //  for DFS
    Region_Node* stack_top;
    Last_point_index = 0;
    // First step using the DFS
    while(!subtree_root.empty()){
        root_stack.push(subtree_root.back());
        subtree_root.pop_back();
    }
    while(root_stack.size() != 0){
        if(root_stack.top()->Partitions.size() ==0 && root_stack.top()->must_save==0){
            stack_top = root_stack.top();
            root_stack.pop();
            if(stack_top->Node_child.size() >= 2){
                root_stack.push(stack_top->Node_child[1]);
                root_stack.push(stack_top->Node_child[0]);
            }
        }
        else{//else if ( root_stack.top() -> Partitions.size() != 0  || root_stack.top()->must_save==1){   // we don't need to travese this subset  treat it as a leaf !!!
            BFS_mark_subroot(root_stack.top());
            subtree_root.push_back(root_stack.top()); // restore the new node for subtree_root
            root_stack.pop();
        }
    }
    for(unsigned i =0; i < Bsp_Tree.size(); i++)
        Bsp_Tree[i]-> must_save=0;
}

void tree::BFS_mark_subroot(Region_Node *& save_region){
    queue <Region_Node *> bfs_queue;

    if(Last_point_index < save_region->Data_point_start){
        for(int i=0; i< save_region->length; i++)
            Point_index[Last_point_index + i] = Point_index[save_region->Data_point_start + i];
        save_region->Data_point_start = Last_point_index;
    }
    Last_point_index += save_region->length;
    bfs_queue.push(save_region);

    while(!bfs_queue.empty()){
        for(unsigned i=2; i<bfs_queue.front()->Node_child.size(); i++)
            bfs_queue.front()->Node_child[i]->must_save=1;
        if(bfs_queue.front()->Node_child.size() >= 2){
            bfs_queue.front()->Node_child[0]->Data_point_start = bfs_queue.front()->Data_point_start;
            bfs_queue.front()->Node_child[1]->Data_point_start = bfs_queue.front()->Data_point_start + bfs_queue.front()->Node_child[0]->length;
            bfs_queue.push(bfs_queue.front()->Node_child[0]);
            bfs_queue.push(bfs_queue.front()->Node_child[1]);
        }
        bfs_queue.pop();
    }
}

void tree::Erase_partition(Region_Node * & region_select, int nIter){

    int suc=0;
    for(unsigned i = 0 ; i < region_select->Partitions.size(); i++){
        if(region_select->Partitions[i]->PartitionsID == Pts[nIter].PartitionsID){
            region_select->Partitions.erase(region_select->Partitions.begin()+i);
            suc=1;
            break;
        }
    }
    if(suc==0){
        cout<<"error not earse successfule in rand "<<endl;
    }
}

void tree::Resample_function(){
#ifdef PROFILE
    if(Dim == 1)
        timer.TimerStart_cpu(timer.sTime_resample_C);
    else
        timer.TimerStart_cpu(timer.sTime_resample_NC);
#endif

    int halfweight = Resampling;
    unsigned size_of_parts= Pts.size();
    vector<double> ec (size_of_parts , 0);
    for(unsigned i = 0 ; i < size_of_parts ; i++) {
        if(halfweight ==2 ){
            Pts[i].w[List_index] /= 2.0;
            ec[i] = exp ( Pts[i].w[List_index] );
        }
        else {
            ec[i]= exp( Pts[i].w[List_index] );
            Pts[i].w[List_index] =1 ;
        }
    }

    int Old_list_index= List_index;
    int New_list_index=  (Old_list_index+1) %2 ; // =old +1 /2   

    for(unsigned i = 0 ; i < size_of_parts ; i ++){
        int r_target= rand_int(ec);
        Pts[i].region[New_list_index]= Pts[r_target].region[Old_list_index];
        Pts[i].cutregs[New_list_index] = Pts[r_target].cutregs[Old_list_index];
        Pts[i].w[New_list_index]=Pts[r_target].w[Old_list_index];
        Pts[i].score[New_list_index] = Pts[r_target].score[Old_list_index];
        if(Dim != 1){
            // removing the Partitions in the region
            for(unsigned j=0; j<Pts[i].region[Old_list_index].size(); j++){
                Erase_partition(Pts[i].region[Old_list_index][j], i);
            }
            for(unsigned j=0; j<Pts[i].region[New_list_index].size(); j++){
                Pts[i].region[New_list_index][j]->Partitions.push_back(&Pts[i]);
            }
        }
    }
    List_index = New_list_index;

#ifdef PROFILE
    if(Dim == 1)
        timer.TimerFinish_cpu(timer.tTime_resample_C , timer.sTime_resample_C);
    else
        timer.TimerFinish_cpu(timer.tTime_resample_NC , timer.sTime_resample_NC);
#endif
}

void tree::SaveTree(){
    cout<<"SaveTree"<<endl;
    /*for(unsigned i=0; i < Bsp_Tree.size(); i++){
        if(Bsp_Tree[i]->must_save != 0){
            cout<<"Bsp_Tree["<<i<<"]!=0"<<endl;
            getchar();
        }
    }*/
    stack <Region_Node *>  root_stack; //  for DFS
    Region_Node* stack_top;
    for(unsigned i=0; i < P_bestlevel.region[List_index_best].size(); i++){
        root_stack.push(P_bestlevel.region[List_index_best][i]);
    }
    while(!root_stack.empty()){
        root_stack.top()->must_save = 1;
        stack_top = root_stack.top();
        root_stack.pop();
        if(stack_top->Node_parent != NULL)
            root_stack.push(stack_top->Node_parent);
    }
    
    root_stack.push(Bsp_Tree[0]);
    while(!root_stack.empty()){
        stack_top = root_stack.top();
        root_stack.pop();
        for(unsigned i=0; i < stack_top->Node_child.size(); i++){
            root_stack.push(stack_top->Node_child[i]); 
        }
        int event=0;
        for(unsigned i=0; i < stack_top->Node_child.size(); i++){
            if(stack_top->Node_child[i]->must_save==1){
                if(stack_top->Node_child[i+1]->must_save!=1){
                    cout<<"Something weird with left child & right child!"<<endl;
                    getchar();
                }
                stack_top->Node_child[0] = stack_top->Node_child[i];
                stack_top->Node_child[1] = stack_top->Node_child[i+1];
                stack_top->BestCut = i/2; 
                event =1;
                break;
            }
        }
        if(event == 1)
            stack_top->Node_child.erase(stack_top->Node_child.begin()+2, stack_top->Node_child.end());
        else
            stack_top->Node_child.clear();
        if(stack_top->must_save == 0)
            stack_top->destroy();
    }
    
    //debug
    queue <Region_Node *> bfs_queue;
    bfs_queue.push(&Root);
    unsigned count=0, count2=0;
    while(!bfs_queue.empty()){
        if(bfs_queue.front()->Node_child.size()!=2){
            //cout<<bfs_queue.front()->Node_child.size()<<endl;
            //getchar();
            count++;
        }
        else
            count2++;
        for(unsigned i=0; i<bfs_queue.front()->Node_child.size(); i++)
            bfs_queue.push(bfs_queue.front()->Node_child[i]);
        bfs_queue.pop();
    }
    if(count!=P_bestlevel.region[List_index_best].size())
        cout<<count<<"!="<<P_bestlevel.region[List_index_best].size()<<endl;
    cout<<count<<" "<<count2<<endl;
    //debug
    /*for(unsigned i=0; i < Bsp_Tree.size(); i++){
        int total=0;
        for(unsigned j=0; j<Bsp_Tree[i]->Node_child.size(); j++){
            total += Bsp_Tree[i]->Node_child[j]->must_save;
        }
        if(total > 2){
            cout<<"Node"<<i<<" has more than 2 must_save child"<<endl;
            getchar();
        }
    }*/

    /*int total=0;
    for(unsigned i=0; i < Bsp_Tree.size(); i++)
        total += Bsp_Tree[i]->must_save;
    cout<<total<<" "<<P_bestlevel.region[List_index_best].size()<<endl;*/
    cout<<"size of Region_Node* = "<<sizeof(Region_Node *)<<endl;
}

double tree::inv_Fprob(double Fx, vector<pair<double,double> >& intervals){
    double x=-100000;
    if(intervals.size()==1) return Fx;
    if( Fx<intervals[0].second)  {
        x=intervals[1].first*Fx/intervals[0].second;
        return x;
    }
    for(unsigned i=1; i<(intervals.size()-1); i++){
        if( Fx<intervals[i].second){
            x = intervals[i].first+(intervals[i+1].first-intervals[i].first)*(Fx-intervals[i-1].second)/(intervals[i].second-intervals[i-1].second);
            return x;
        }
    }
    unsigned i= intervals.size()-1;
    x = intervals[i].first+(1-intervals[i].first)*(Fx-intervals[i-1].second)/(1-intervals[i-1].second);
    return x;

}


double tree::inv_Ftransform_for_output(double Fx ,     vector <pair<double , double> > intervals   ){

    double x = inv_Fprob(Fx , intervals);
    return x;
}

vector<double> tree::inv_Ftransform_t(vector<double> Fx,   Partition &P    /*int bestP */  ){
    vector<double> x(Fx.size(), 0);
    vector<pair<double,double> > intervals = onedimpartitiontointerval_t_copy( P);
    for(unsigned i=0; i< Fx.size(); i++)  x[i] = inv_Fprob(Fx[i], intervals);
    return x;
}

vector <double> tree::Ftransform_t(vector<double> x ,  Partition &P /* int bestP*/   ) {

    vector <double>F(x.size() , 0 );
    vector <pair<double,double> > intervals =  onedimpartitiontointerval_t_copy(P);
    for(unsigned i=0 ; i < x.size(); ++i)  
    {
        if(i<10){
            cerr<<"Intervals " <<intervals[0].first
                <<"," <<intervals[0].second <<'\n';
        }

        F[i]=Fprob(x[i] , intervals);
    }
    return F;
}

// transform the partition to the pair<double,double> recording the starting point and the cdf. order the intervals by the position. the first pair.first=0
vector<pair<double,double> > tree::onedimpartitiontointerval_t_copy  (Partition  &P ){
    unsigned numreg = P.region[List_index].size();

    pair<double,double> range;
    int totalnum=0;
    for(unsigned i=0; i< numreg; i++){totalnum += P.region[List_index][i]->length;  }  // +

    vector<pair<double,double> > intervals;
    if(numreg==1) {
        range.first=0;
        range.second=1;
        intervals.push_back(range);
        return (intervals);
    }
    for(unsigned i=0; i< numreg; i++){
        range.first = convert_ranges(P.region[List_index][i]->Region_mask[0]).first;   // Andy
        range.second = (double)(P.region[List_index][i]->length)/totalnum/exp(lprod_usint_mask(P.region[List_index][i]->Region_mask ));                   // Andy
        intervals.push_back(range);
    }

    sort( intervals.begin(), intervals.end() , compareinterval);


    //  transform the intervals[i].second to the cdf.
    intervals[0].second =  (intervals[1].first- intervals[0].first)* intervals[0].second;
    for(unsigned i=1; i<(numreg-1); i++){
        intervals[i].second = intervals[i-1].second + (intervals[i+1].first- intervals[i].first)* intervals[i].second;
    }
    intervals[numreg-1].second = intervals[numreg-2].second + (1- intervals[numreg-1].first)* intervals[numreg-1].second;

    if(fabs(intervals[numreg-1].second -1) >0.00001){
        cout<<"wrong in onedimpartitiontointerval\n";
        for(unsigned i=0; i<numreg; i++){
            cout<<intervals[i].first<<" "<<intervals[i].second<<"\n";
        }
        exit(9);
    }
    return intervals;
}

double tree::Fprob(double x, vector<pair<double,double> >& intervals){
    double y=0;
    bool cont=true;
    if(intervals.size()==1) return x;
    if( x<intervals[1].first)  {
        y=intervals[0].second*x/intervals[1].first;
        cont=false;
    }
    for(unsigned i=1; cont && i<(intervals.size()-1); i++){
        if( x<intervals[i+1].first){
            y = intervals[i-1].second+(intervals[i].second-intervals[i-1].second)*(x-intervals[i].first)/(intervals[i+1].first-intervals[i].first);
            cont=false;
        }
    }
    if(cont) {
        unsigned i= intervals.size()-1;
        y = intervals[i-1].second+(1-intervals[i-1].second)*(x-intervals[i].first)/(1-intervals[i].first);
    }
    if(y>=-0.0001 && y<0.0000001)  y=0.000001;         // make sure 0 and 1 don't appear.
    else if(y>0.999999 && y <= 1.0001)  y=0.999999;
    else if (y<-0.0001 || y>1.0001){
        cerr << "y: " << y << '\n';
        cerr << "wrong y in Fprob!";
        exit(10);
    }
    return y;

}

/*
bool compareinterval (const pair<double, double> &i, const pair<double, double> &j) {
    return (i.first<j.first);
}
*/

// transform the partition to the pair<double,double> recording the starting point and the cdf. order the intervals by the position. the first pair.first=0
vector<pair<double,double> > tree::onedimpartitiontointerval_t( int bestP ){
    unsigned numreg = Pts[bestP].region[List_index].size(); 
    pair<double,double> range;
    int totalnum=0;
    for(unsigned i=0; i< numreg; i++){totalnum += Pts[bestP].region[List_index][i]->length;  }  // +
    vector<pair<double,double> > intervals;
    if(numreg==1) {
        range.first=0;
        range.second=1;
        intervals.push_back(range);
        return (intervals);
    }
    for(unsigned i=0; i< numreg; i++){
        range.first = convert_ranges(Pts[bestP].region[List_index][i]->Region_mask[0]).first;   // Andy 
        range.second = (double)(Pts[bestP].region[List_index][i]->length)/totalnum/exp(lprod_usint_mask(Pts[bestP].region[List_index][i]->Region_mask ));                   // Andy 
        intervals.push_back(range);
    }
    sort( intervals.begin(), intervals.end() , compareinterval);
    //  transform the intervals[i].second to the cdf.
    intervals[0].second =  (intervals[1].first- intervals[0].first)* intervals[0].second;
    for(unsigned i=1; i<(numreg-1); i++){
        intervals[i].second = intervals[i-1].second + (intervals[i+1].first- intervals[i].first)* intervals[i].second;
    }
    intervals[numreg-1].second = intervals[numreg-2].second + (1- intervals[numreg-1].first)* intervals[numreg-1].second;

    if(fabs(intervals[numreg-1].second -1) >0.00001){
        cout<<"wrong in onedimpartitiontointerval\n";
        for(unsigned i=0; i<numreg; i++){
            cout<<intervals[i].first<<" "<<intervals[i].second<<"\n";
        }
        exit(9);
    }
    return intervals;
}

// from usint_mask to a pair<double, double> range with its small and large boundary
pair<double, double> tree::convert_ranges(usint_mask root) {
    double up = 0.5;
    //          double low=0;
    unsigned int u = root.x;
    unsigned int v = root.mask;
    while (v > 0) {
        bool ind = (u & 1u);
        if (ind) {
            up += 0.5;
            //                          low += 0.5;
        }
        //                      low/=2;
        up /= 2;
        u = u >> 1;
        v = v >> 1;
    }
    up = up * 2;
    double low = up - 1.0 / (double) (root.mask + 1);
    pair<double, double> range;
    range.first = low;
    range.second = up;
    return range;
}

double tree::lprod_usint_mask(const vector<usint_mask>& reg_code) {                 // return log of the area
    double result = 0;
    for (unsigned i = 0; i < reg_code.size(); i++) {
        result -= log((double)(reg_code[i].mask+1));
    }
    return result;
}

void tree::print_partition1_t(const Partition &P) {
    int totalnum = 0;

    for (unsigned i = 0; i < P.region[List_index].size(); i++) {
        totalnum += P.region[List_index][i]->length;
    }

    for (unsigned i = 0; i < P.region[List_index].size(); i++) {
        pair<double, double> range;
        for (unsigned d = 0; d < P.region[List_index][0]->Region_mask.size(); d++){ // ??? why it is reg 0 ; 
            range = convert_ranges(P.region[List_index][i]->Region_mask[d]);
            cout << range.first << " " << range.second << " ";
        }
        cout <<  (double) P.region[List_index][i]->length / totalnum / exp(lprod_usint_mask(P.region[List_index][i]->Region_mask)) << " " << P.region[List_index][i]->length << endl;
    }

}

void tree::print_partition2_t( ofstream  &outfile    ,  const Partition &P) {
    int totalnum = 0;

    for (unsigned i = 0; i < P.region[List_index].size(); i++) {
        totalnum += P.region[List_index][i]->length;
    }

    for (unsigned i = 0; i < P.region[List_index].size(); i++) {
        pair<double, double> range;
        for (unsigned d = 0; d < P.region[List_index][0]->Region_mask.size(); d++){ // ??? why it is reg 0 ;
            range = convert_ranges(P.region[List_index][i]->Region_mask[d]);
            outfile << setprecision(12) << range.first << " " << setprecision(12) << range.second << " ";
        }
        outfile << setprecision(12) << (double) P.region[List_index][i]->length / totalnum / exp(lprod_usint_mask(P.region[List_index][i]->Region_mask)) << " " << P.region[List_index][i]->length << endl;
    }

}

void tree::print_partition3_t(ofstream &outfile, const Partition &P, const vector<double> & mmax, const vector<double>& mmin) {
    int totalnum = 0;

    for (unsigned i = 0; i < P.region[List_index_best].size(); i++) {
        totalnum += P.region[List_index_best][i]->length;
    }

    for (unsigned i = 0; i < P.region[List_index_best].size(); i++) {
        pair<double, double> range;
        double newarea=1.0;
        for (unsigned d = 0; d < P.region[List_index_best][0]->Region_mask.size(); d++){
            range = convert_ranges(P.region[List_index_best][i]->Region_mask[d]);
            double newrange1 = range.first *(mmax[d] - mmin[d])+ mmin[d];
            double newrange2 = range.second *(mmax[d] - mmin[d])+ mmin[d];
            newarea = newarea*(newrange2-newrange1);
            outfile<<newrange1<<" "<<newrange2<< " ";
        }
        outfile<< (double) P.region[List_index_best][i]->length / (double)totalnum/newarea<<" "<<P.region[List_index_best][i]->length<<endl;
    }
}

void tree::output_tree(ofstream &outfile, bool normalization){
    stack <Region_Node *> reg_stack;
    pair<double,double> reg_range;
    Region_Node* stack_top;
    double amount = 0;
    int node_id = 0;
    int leaf_count = 0;
    double J = 1;
    
    Root.NodeID = node_id;
    node_id++;
    Root.volume = Root.find_volume();
    Root.prob = (double)Root.length/(double)Samplesize;
    Root.cut_dimension = Root.cut_dim();
    
    for (unsigned int i=0; i < Root.Node_child.size(); i++){        
        reg_stack.push(Root.Node_child[i]);    
    }
    
    while (!reg_stack.empty()){        
        stack_top = reg_stack.top();        
        stack_top->NodeID = node_id;        
        stack_top->volume = stack_top->find_volume();
        stack_top->prob = (double)stack_top->length/(double)Samplesize;                
        stack_top->cut_dimension = stack_top->cut_dim();                
        node_id++;                        
        reg_stack.pop();        
        for (unsigned int i=0; i < stack_top->Node_child.size(); i++){            
            reg_stack.push(stack_top->Node_child[i]);       
        }        
        if (stack_top->Node_child.size() == 0){	//Leaf found            
            leaf_count++;            
            amount += stack_top->length;        
        }    
    }
    // Output //
    outfile << Root.NodeID <<'\t';
    for (int d=0; d<Dim; d++){
        reg_range = convert_ranges(Root.Region_mask[d]);
        if (normalization){
            reg_range.first = mmin[d] + (mmax[d] - mmin[d])*reg_range.first;
            reg_range.second = mmin[d] + (mmax[d] - mmin[d])*reg_range.second;
            J *= (mmax[d] - mmin[d]);
        }
        //outfile << setprecision(12) << reg_range.first <<'\t'<< setprecision(12) << reg_range.second <<'\t';
        outfile << reg_range.first <<'\t'<< reg_range.second <<'\t';
    }
    outfile << Root.prob/Root.volume/J <<'\t';
    outfile << Root.length <<'\t';    
    outfile << Root.cut_dimension<<'\t';    
    outfile << -1 <<'\t';    
    outfile << Root.Node_child[0]->NodeID <<'\t';    
    outfile << Root.Node_child[1]->NodeID <<endl;
    
    for (unsigned int i=0; i < Root.Node_child.size(); i++){        
        reg_stack.push(Root.Node_child[i]);            
    }
    while (!reg_stack.empty()){        
        stack_top = reg_stack.top();        
        outfile << stack_top->NodeID <<'\t';        
        for (int d=0; d<Dim; d++){
            reg_range = convert_ranges(stack_top->Region_mask[d]);
            if (normalization){
                reg_range.first = mmin[d] + (mmax[d] - mmin[d])*reg_range.first;
                reg_range.second = mmin[d] + (mmax[d] - mmin[d])*reg_range.second;
            }
            //outfile << setprecision(12) << reg_range.first <<'\t'<< setprecision(12) << reg_range.second<<'\t';
            outfile << reg_range.first <<'\t'<< reg_range.second <<'\t';
        }
        outfile << stack_top->prob/stack_top->volume/J <<'\t';
        outfile << stack_top->length <<'\t';        
        outfile << stack_top->cut_dimension<<'\t';        
        outfile << stack_top->Node_parent->NodeID <<'\t';
        
        reg_stack.pop();
        for (unsigned int i=0; i < stack_top->Node_child.size(); i++){            
            reg_stack.push(stack_top->Node_child[i]);        
        }
        if (stack_top->Node_child.size()==0){            
            outfile << -1 <<'\t';            
            outfile << -1 <<endl;        
        }        
        else{            
            outfile << stack_top->Node_child[0]->NodeID <<'\t';            
            outfile << stack_top->Node_child[1]->NodeID <<endl;        
        }
    }
}