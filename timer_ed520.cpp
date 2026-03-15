#include <iostream>
#include "timer_ed520.h"
using namespace std;

void TIMER::Begin(){
    tProgramBegin = clock();
}

void TIMER::End(){
    //return (clock() - tProgramBegin)/CLOCKS_PER_SEC;
    //return (clock() - tProgramBegin);
    ExecutionTime = (clock() - tProgramBegin);
}


void TIMER::TimerStart_cpu( float& sTime  ){    // using the clock timer to measure the time

     sTime = clock();
}

void TIMER::TimerFinish_cpu( float & tTime , float & sTime   ){    // using the colck timer to measure the time

     tTime=tTime + (clock()- sTime );
}



/*void TIMER::TimerStart(){     // using the CUDA time function (very percise )
    cudaEventRecord(start, 0);
    //start1 = clock();
}

void TIMER::TimerFinish(float& tTime){
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    tTime = tTime + elapsedTime;
    //tTime = tTime + (clock() - start1);
}*/

float TIMER::RunTimeProfile(){
    cout<<"====== Runtime (ms) ======"<<endl;
    cout<<"Total Execution Time = "<<ExecutionTime/CLOCKS_PER_SEC<<" seconds"<<endl;
#ifdef PROFILE
    cout<<"C_SIS time is :\t\t\t" <<tTime_C_SIS*1000/CLOCKS_PER_SEC<<endl;
    cout<<"NC_SIS time is :\t\t"<<tTime_NC_SIS*1000/CLOCKS_PER_SEC<<endl;
    cout<<"init_C time is :\t\t"<<tTime_init_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"resample_C time is :\t\t"<<tTime_resample_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"get_info_C time is :\t\t"<<tTime_get_info_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"rand_choose_C time is :\t\t"<<tTime_rand_choose_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"Count_C time is :\t\t"<<tTime_Count_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"Calculate_weight_C time is :\t"<<tTime_Calculate_weight_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"Two_Maxc_in_one_part_C time is :"<<tTime_Two_Maxc_in_one_part_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"init_NC time is :\t\t"<<tTime_init_NC*1000/CLOCKS_PER_SEC<<endl;
    cout<<"resample_NC time is :\t\t"<<tTime_resample_NC*1000/CLOCKS_PER_SEC<<endl;
    cout<<"get_info_NC time is :\t\t"<<tTime_get_info_NC*1000/CLOCKS_PER_SEC<<endl;
    cout<<"rand_choose_NC time is :\t"<<tTime_rand_choose_NC*1000/CLOCKS_PER_SEC<<endl;
    cout<<"Count_NC time is :\t\t"<<tTime_Count_NC_cpu*1000/CLOCKS_PER_SEC<<endl;
    cout<<"Calculate_weight_NC time is :\t"<<tTime_Calculate_weight_NC*1000/CLOCKS_PER_SEC<<endl;
    cout<<"Two_Maxc_in_one_part_NC time is :"<<tTime_Two_Maxc_in_one_part_NC*1000/CLOCKS_PER_SEC<<endl;
    cout<<"Count Max BPscore is :\t\t"<<tTime_BPscore*1000/CLOCKS_PER_SEC<<endl;
    /*cout<<"Cuda initilization_C is :\t\t" <<tTime_cuda_init_C*1000/CLOCKS_PER_SEC<<endl;
    cout<<"CUDA_CreateDataArray_Dynamic_C time:\t"<<tCUDA_CreateDataArray_Dynamic_C<<endl;
    cout<<"MemoryCopy_H2D_Dynamic time_C:\t\t"<<tMemoryCopy_H2D_Dynamic_C<<endl;
    cout<<"Count_GPU_C time:\t\t\t"<<tCount_GPU_C<<endl;
    cout<<"Reduction_C time:\t\t\t"<<tReduction_C<<endl;
    cout<<"MemoryCopy_D2H_C time:\t\t\t"<<tMemoryCopy_D2H_C<<endl;
    
    cout<<"Cuda initilization_NC is :\t\t" <<tTime_cuda_init_NC*1000/CLOCKS_PER_SEC<<endl;
    cout<<"CUDA_CreateDataArray_Dynamic_NC time:\t"<<tCUDA_CreateDataArray_Dynamic_NC<<endl;
    cout<<"MemoryCopy_H2D_Dynamic_NC time:\t\t"<<tMemoryCopy_H2D_Dynamic_NC<<endl;
    cout<<"Count_GPU_NC time:\t\t\t"<<tCount_GPU_NC<<endl;
    cout<<"Reduction_NC time:\t\t\t"<<tReduction_NC<<endl;
    cout<<"MemoryCopy_D2H_NC time:\t\t\t"<<tMemoryCopy_D2H_NC<<endl;
    
    cout<<"number of calling the cuda:\t\t" <<number_in_cuda<<endl;*/
#endif
	return ExecutionTime/CLOCKS_PER_SEC;
}
