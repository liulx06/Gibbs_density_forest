#ifndef TIMER_ED520_H
#define TIMER_ED520_H

class TIMER{
    public:
    float tProgramBegin;
    float ExecutionTime;
    void Begin();
    void End();
    TIMER(){
        tProgramBegin = 0;
        ExecutionTime = 0;
        /// CUDA /// 
        /*cudaEventCreate (&start);
        cudaEventCreate (&stop);
        tCUDA_CreateDataArray_Dynamic_NC=0;
        tCUDA_CreateDataArray_Dynamic_C=0;
        tMemoryCopy_H2D_Dynamic_NC=0;
        tMemoryCopy_H2D_Dynamic_C=0;
        tCount_GPU_NC=0;
        tCount_GPU_C=0;
        tReduction_NC=0;
        tReduction_C=0;
        tMemoryCopy_D2H_NC=0;
        tMemoryCopy_D2H_C=0;
        number_in_cuda=0;*/
        /// CPU ///
        sTime_cuda_init_C=0; tTime_cuda_init_C=0;
        sTime_cuda_init_NC=0; tTime_cuda_init_NC=0;
        sTime_C_SIS=0; tTime_C_SIS=0;
        sTime_NC_SIS=0; tTime_NC_SIS=0;
        sTime_resample_C=0; tTime_resample_C=0;
        sTime_resample_NC=0; tTime_resample_NC=0;
        sTime_get_info_C=0; tTime_get_info_C=0;
        sTime_get_info_NC=0; tTime_get_info_NC=0;
        sTime_rand_choose_C=0; tTime_rand_choose_C=0;
        sTime_rand_choose_NC=0; tTime_rand_choose_NC=0;
        //add by CY
        sTime_Count_C=0;tTime_Count_C=0;
        sTime_Count_NC_cpu=0;tTime_Count_NC_cpu=0;
        sTime_Count_NC_gpu=0;tTime_Count_NC_gpu=0;
        sTime_Calculate_weight_C=0;tTime_Calculate_weight_C=0;
        sTime_Calculate_weight_NC=0;tTime_Calculate_weight_NC=0;
        sTime_Two_Maxc_in_one_part_C=0;tTime_Two_Maxc_in_one_part_C=0;
        sTime_Two_Maxc_in_one_part_NC=0;tTime_Two_Maxc_in_one_part_NC=0;
        sTime_init_C=0;tTime_init_C=0;
        sTime_init_NC=0;tTime_init_NC=0;
        sTime_CUDA_INIT=0;tTime_CUDA_INIT=0;
        sTime_BPscore=0;tTime_BPscore=0;
        CPU_COUNT=0; GPU_COUNT=0; HIT_COUNT=0;
    }
    float RunTimeProfile();
    ////// define Profile2 ///////
    /*cudaEvent_t start;
    cudaEvent_t stop;
    float tCUDA_CreateDataArray_Dynamic_NC;
    float tCUDA_CreateDataArray_Dynamic_C;
    float tMemoryCopy_H2D_Dynamic_NC;
    float tMemoryCopy_H2D_Dynamic_C;
    float tCount_GPU_NC;
    float tCount_GPU_C;
    float tReduction_NC;
    float tReduction_C;
    float tMemoryCopy_D2H_NC;
    float tMemoryCopy_D2H_C;
    void TimerStart();
    void TimerFinish(float& tTime);

    int number_in_cuda;*/

    float sTime_cuda_init_C, tTime_cuda_init_C;
    float sTime_cuda_init_NC, tTime_cuda_init_NC;
    float sTime_C_SIS, tTime_C_SIS;
    float sTime_NC_SIS, tTime_NC_SIS;
    float sTime_resample_C, tTime_resample_C;
    float sTime_resample_NC, tTime_resample_NC;
    float sTime_get_info_C, tTime_get_info_C;
    float sTime_get_info_NC, tTime_get_info_NC;
    float sTime_rand_choose_C, tTime_rand_choose_C;
    float sTime_rand_choose_NC, tTime_rand_choose_NC;
    //add by CY
    float sTime_Count_C,tTime_Count_C;
    float sTime_Count_NC_cpu,tTime_Count_NC_cpu;
    float sTime_Count_NC_gpu,tTime_Count_NC_gpu;
    float sTime_Calculate_weight_C,tTime_Calculate_weight_C;
    float sTime_Calculate_weight_NC,tTime_Calculate_weight_NC;
    float sTime_Two_Maxc_in_one_part_C,tTime_Two_Maxc_in_one_part_C;
    float sTime_Two_Maxc_in_one_part_NC,tTime_Two_Maxc_in_one_part_NC;
    float sTime_init_C,tTime_init_C;
    float sTime_init_NC,tTime_init_NC;
    float sTime_CUDA_INIT,tTime_CUDA_INIT;
    float sTime_BPscore,tTime_BPscore;
    void TimerStart_cpu(float& sTime);
    void TimerFinish_cpu(float& tTime, float& sTime);
    //random choose profile
    float sTime_rand_2,tTime_rand_2;
    float sTime_rand_3,tTime_rand_3;
    int CPU_COUNT;
    int GPU_COUNT;
    int HIT_COUNT;
}; 
#endif

