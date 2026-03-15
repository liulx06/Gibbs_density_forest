#ifndef BSP_Tree_H
#define BSP_Tree_H
#include "stl.h"
#include "header.h"

#define SCALE 1

bool compareinterval (const pair<double, double> &i, const pair<double, double> &j);

class Region_Node;

class Partition{
    public:
        Partition();
        void destroy();
        vector <Region_Node *> region_1;    // it is one region nood chain  
        vector <Region_Node *> region_2;    // it is another region node chain
        vector < vector< Region_Node *> > region;  // ??? it may be the tree structure for fast cut searching 
        
        Region_Node* region_replace;       //Replaced  region at last cut
        double score[2];    //BPscore
        double w[2];
        int cutregs[2]; // 0 for the one chain, 1 for another chain 
        int first;      // we use the size of the region[List_index] == 1 to replace the first flag  
        double cmax ; 
        double cmax2;
        int PartitionsID;
            
};
class Region_Node{
    public:
        Region_Node(){
            must_save=0;
        }
        int must_save;
        vector<Partition*> Partitions;  // it would be important later 
        // boundary   
        usint_mask um;
        vector<usint_mask> Region_mask;
        // end of boundary 
        bool * unique_dim;
        vector<double> c_original;
        vector<double> c_modified;
        vector<double> c_local_max;        //record the maximum c_original & 2nd maximum c_original
        //double c_local_max;
        int Data_point_start;           //start index of this region 
        // modified int to double
        double length;                  //modfied 412: # of data points of this region
        double prob;			//+ 412 for CDF operation usage
        double volume;			//+ 412 for CDF operation usage
        //
        vector <int> sub_length;
        vector <int> Dim_has_choose; 
        //  we need it to kown that whether we are going to modify the index _amount ! use it in the rand_choose func
        vector< Region_Node * >Node_child;//ignore the left and right child reference is suitable !!! 
        Region_Node* Node_parent;
        int BestCut;
        int NodeID;		//+ 412
        int cut_dimension;	//+ 412
        int cut_dim();		//+ 412
        double find_volume(); 	//+ 412
        void destroy();
};
class tree{
    public:     
        int Dim;
        double *Data_array;
        int  *Point_index;        //  samplenum * partition 
        
        int Last_point_index;    //  to record the 
        vector<vector<double> > Data;// Lu's we can replace this one later ! 
        double maxpercentage;
        int n;                   //  total amount partitions 
        int Levels;              //  total cut amount 
        int Resampling;          //  resample enabled or not 
        int Steps;               //  resample per step  
        int Samplesize;
        Region_Node Root;        //  root of BSP tree
        vector<Region_Node *> Bsp_Tree;  // Root Region it contain all the region 
        vector<Region_Node *> subtree_root;   // for Data_point_index rearrange, 
        // the following are used in the SIS_tree func
        int BestPartID;      
        vector<Partition> Pts;   //  total partition  
        bool discrete;           //  constructer assign the value
        double beta;             //  constructer assign the value
        Partition P_bestlevel;
        vector<double> mmax;	//+ 412
        vector<double> mmin;	//+ 412
        
        void tree_para_init( );
        void tree_para_destroy();
        void tree_nonpara_init( );
        void para_share_init();


        // SIS function and the variable that the function need
        void SIS_tree();   

        bool is_milestone(int i , int plevels );

        int modelBP;   
        double currlBPscore; 
        double maxlBPscore;
        int maxlevel;
        bool exhaust;

        int List_index;
        int List_index_best;     //  record the best partition of its List_index
        void Resample_function();

        //   sample_tree and it's necessary variable 
        bool Sample_tree();
        double smooth;
    
        void get_info_from_onetree(int  nIter );
        void Two_Maxc_in_one_part( Partition & P ); // find the max ans second in a partition 

        int *Tmp_count;
        bool *Tmp_unique;
        
        int hit_in_count;   // a flag for calculate  
        void Count( Region_Node * region_cut  );
        void Calculate_weight(int, int); 
        void rand_choose_one_nextleveltree(int nIter);
        void Region_Create( Region_Node *region_cut, Region_Node *&New_Region_L, Region_Node *&New_Region_R, int i);
        void Erase_partition( Region_Node * & region_select, int nIter);
        void Data_point_defragment();  // for Data_point_rearrangement
        void BFS_mark_subroot(Region_Node* & save_region);
    
        void SaveTree();

        //for statstic function 
        double getmaxlBP_t(vector<Partition>& newPs);
        double lBPscore_t(Partition &P);  // it would be called in the getmaxlBP 
        // for transform to better 
        vector <double> Ftransform_t( vector<double> x , Partition &P   /*  int bestP */ );
        double Fprob( double x , vector<pair<double,double> > & intervals);
        vector<pair<double,double> > onedimpartitiontointerval_t(int bestP);
        vector<pair<double,double> > onedimpartitiontointerval_t_copy(Partition &P);

        bool comparaintterval( const pair<double , double> &i ,const pair<double,double> &j);
        pair <double ,double> convert_ranges (usint_mask root );
        double lprod_usint_mask( const vector<usint_mask> & reg_code );
        double inv_Fprob(double Fx , vector<pair<double ,double> > &intervals);
        vector<double> inv_Ftransform_t( vector<double> Fx ,   Partition &P    /* int bestP */ );
        double inv_Ftransform_for_output(double Fx,vector <pair<double , double> > intervals /*  Partition &P  */ );
        // for print
        void print_partition1_t(const Partition &P);
        void print_partition2_t(ofstream& outfile, const Partition &P);
        void print_partition3_t(ofstream& outfile, const Partition &P, const vector<double>& mmax, const vector<double>& mmin);
        void output_tree(ofstream& outfile, bool normalization);
};

#endif
