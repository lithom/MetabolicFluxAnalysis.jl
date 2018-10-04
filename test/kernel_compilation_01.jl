

using OpenCL

str_kernel_test_01 = "

__kernel void fwsim(__global const float *semus,
                     __global const float *fluxes,
                     __global float *mds_out)
   {
     const int gid = get_group_id(0);
     const int idx = get_local_id(0);
     const int idy = get_local_id(1);

     //int offset_semus   = 0;
     //int offset_fluxes = gid * 17;
     //int offset_out    = gid * 3;

     // Here we store the sizes of the given matrices in a step
     int mAm = 0;
     int mAn = 0;
     int mBm = 0;
     int mBn = 0;
     int mXm = 0;
     int mXn = 0;
     int mYm = 0;
     int mYn = 0;




__local float mA[25];
__local float mA_L[25];
__local float mB[10];
__local float mX[12];
__local float mY[8];
__local float mBY[12];


//Required EMRs: 11  Required Fluxes: 6  Target MDs: 4
int offset_emrs   = 0;
int offset_fluxes = gid * 6;
int offset_target = gid * 4;
//-----Init EMRs-------------------------------------------------------------------------------------
float mr_sA_2_M_0 = semus[offset_emrs + 1];
float mr_sA_2_M_1 = semus[offset_emrs + 2];
float mr_sA_3_M_0 = semus[offset_emrs + 3];
float mr_sA_3_M_1 = semus[offset_emrs + 4];
float mr_sA_2_3_M_0 = semus[offset_emrs + 5];
float mr_sA_2_3_M_1 = semus[offset_emrs + 6];
float mr_sA_2_3_M_2 = semus[offset_emrs + 7];
float mr_sA_1_2_3_M_0 = semus[offset_emrs + 8];
float mr_sA_1_2_3_M_1 = semus[offset_emrs + 9];
float mr_sA_1_2_3_M_2 = semus[offset_emrs + 10];
float mr_sA_1_2_3_M_3 = semus[offset_emrs + 11];
//-----Init Fluxes-----------------------------------------------------------------------------------
float v_R01 = fluxes[ offset_fluxes + 1];
float v_R02 = fluxes[ offset_fluxes + 2];
float v_R03 = fluxes[ offset_fluxes + 3];
float v_R04 = fluxes[ offset_fluxes + 4];
float v_R05 = fluxes[ offset_fluxes + 5];
float v_R06 = fluxes[ offset_fluxes + 6];
//---------------------------------------------------------------------------------------------------
//------Declare All Temporary EMR Float values--------------------------------------------------------------
float mr_sC_1_M_1;
float mr_sD_3_M_0;
float mr_sB_1_2_3_M_1;
float mr_sD_2_M_1;
float mr_sF_1_2_3_M_1;
float mr_sD_2_3_M_0;
float mr_sD_1_2_3_M_0;
float mr_sF_1_2_3_M_2;
float mr_sB_3_M_1;
float mr_sB_2_3_M_0;
float mr_sB_1_2_3_M_2;
float mr_sB_2_M_0;
float mr_sD_1_2_3_M_3;
float mr_sB_2_M_1;
float mr_sD_1_2_3_M_1;
float mr_sD_2_3_M_1;
float mr_sB_3_M_0;
float mr_sB_2_3_M_1;
float mr_sF_1_2_3_M_0;
float mr_sD_2_M_0;
float mr_sB_1_2_3_M_3;
float mr_sB_2_3_M_2;
float mr_sD_2_3_M_2;
float mr_sD_1_2_3_M_2;
float mr_sC_1_M_0;
float mr_sF_1_2_3_M_3;
float mr_sD_3_M_1;
float mr_sB_1_2_3_M_0;
//---------------------------------------------------------------------------------------------------
//----Init A---------------------------------------------------------------------------------
if(idx==0 && idy==0){
mAm = 5;
mAn = 5;
mA[0] = ( (  -1.0f  * v_R01   )  +(  -1.0f  * v_R03   )   ) ;
mA[1] = 0.0f;
mA[2] =  v_R04 ;
mA[3] =  v_R02 ;
mA[4] = 0.0f;
mA[5] = 0.0f;
mA[6] = ( (  -1.0f  * v_R01   )  +(  -1.0f  * v_R03   )   ) ;
mA[7] = 0.0f;
mA[8] =  v_R05 ;
mA[9] =  v_R02 ;
mA[10] = 0.0f;
mA[11] = 0.0f;
mA[12] = ( (  -1.0f  * v_R04   )   ) ;
mA[13] = 0.0f;
mA[14] =  v_R05 ;
mA[15] =  v_R03 ;
mA[16] = 0.0f;
mA[17] = 0.0f;
mA[18] = ( (  -1.0f  * v_R02   )  +(  -1.0f  * v_R05   )   ) ;
mA[19] = 0.0f;
mA[20] = 0.0f;
mA[21] =  v_R03 ;
mA[22] = 0.0f;
mA[23] = 0.0f;
mA[24] = ( (  -1.0f  * v_R05   )  +(  -1.0f  * v_R02   )   ) ;
}
//----Init B---------------------------------------------------------------------------
if(idx==1 && idy==0){
mBm = 5;
mBn = 2;
mB[0] = 0.0f;
mB[1] = (  -1.0f  * v_R01   ) ;
mB[2] = 0.0f;
mB[3] = 0.0f;
mB[4] = 0.0f;
mB[5] = (  -1.0f  * v_R01   ) ;
mB[6] = 0.0f;
mB[7] = 0.0f;
mB[8] = 0.0f;
mB[9] = 0.0f;
}
//----Init Y---------------------------------------------------------------------------
if(idx==2 && idy==0){
mYm = 2;
mYn = 2;
mY[0] = ( (  mr_sA_3_M_0   )   ) ;
mY[1] = ( (  mr_sA_2_M_0   )   ) ;
mY[2] = ( (  mr_sA_3_M_1   )   ) ;
mY[3] = ( (  mr_sA_2_M_1   )   ) ;
}

for(int i=0;i<mAm*mAn;i++){ printf(\" %f \",mA[i]); }
printf(\" \\n \")
for(int i=0;i<mBm*mBn;i++){ printf(\" %f \",mB[i]); }
printf(\" \\n \")
for(int i=0;i<mYm*mYn;i++){ printf(\" %f \",mY[i]); }
printf(\" \\n \")

//------COMPUTE!-----------------------------------------------------------------------
//LU-DECOMPOSITION (in place): mA_L * mA = mA
if( idx<5 && idy<5 ){ mA_L[ idx+5*idy ] = (idx==idy)?1.0f:0.0f; }
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 0 / 4
if( idx==0 && idy > 0 && idy < 5 ){
mA_L[ idy + 0 ] = mA[ idy + 0 ] / mA[ 0 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>0 && idx < 5 && idy < 5 ){
mA[ idx + 5*idy ] -= mA_L[ idx + 0 ] * mA[0+ idy*5 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 1 / 4
if( idx==0 && idy > 1 && idy < 5 ){
mA_L[ idy + 5 ] = mA[ idy + 5 ] / mA[ 6 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>1 && idx < 5 && idy < 5 ){
mA[ idx + 5*idy ] -= mA_L[ idx + 5 ] * mA[1+ idy*5 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 2 / 4
if( idx==0 && idy > 2 && idy < 5 ){
mA_L[ idy + 10 ] = mA[ idy + 10 ] / mA[ 12 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>2 && idx < 5 && idy < 5 ){
mA[ idx + 5*idy ] -= mA_L[ idx + 10 ] * mA[2+ idy*5 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 3 / 4
if( idx==0 && idy > 3 && idy < 5 ){
mA_L[ idy + 15 ] = mA[ idy + 15 ] / mA[ 18 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>3 && idx < 5 && idy < 5 ){
mA[ idx + 5*idy ] -= mA_L[ idx + 15 ] * mA[3+ idy*5 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 4 / 4
if( idx==0 && idy > 4 && idy < 5 ){
mA_L[ idy + 20 ] = mA[ idy + 20 ] / mA[ 24 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>4 && idx < 5 && idy < 5 ){
mA[ idx + 5*idy ] -= mA_L[ idx + 20 ] * mA[4+ idy*5 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION DONE!
//----------------------------------
//MATRIX MUL: mBY = mB * mY
//Work item (idx,idy) compute mBY[idx,idy]
if( idx < 5 && idy < 2 ) {
float mulsum = 0.0f;
mulsum += mB[ 5*0 + idx ] * mY[ idy * 2 + 0 ];
mulsum += mB[ 5*1 + idx ] * mY[ idy * 2 + 1 ];
mBY[idx+5*idy] = mulsum;
}
//MATRIX MUL COMPLETE!
barrier(CLK_LOCAL_MEM_FENCE);
// Solve L System
// Solve L System - Row 0
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
mBY[ 0 + 5*idy ] = (mBY[ 0 + 5*idy ] - sum_i ) / mA_L[ 0 ] ;
}
// Solve L System - Row 1
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 0 + 5*idy ];
mBY[ 1 + 5*idy ] = (mBY[ 1 + 5*idy ] - sum_i ) / mA_L[ 6 ] ;
}
// Solve L System - Row 2
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 0 + 5*idy ];
sum_i += mBY[ 1 + 5*idy ];
mBY[ 2 + 5*idy ] = (mBY[ 2 + 5*idy ] - sum_i ) / mA_L[ 12 ] ;
}
// Solve L System - Row 3
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 0 + 5*idy ];
sum_i += mBY[ 1 + 5*idy ];
sum_i += mBY[ 2 + 5*idy ];
mBY[ 3 + 5*idy ] = (mBY[ 3 + 5*idy ] - sum_i ) / mA_L[ 18 ] ;
}
// Solve L System - Row 4
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 0 + 5*idy ];
sum_i += mBY[ 1 + 5*idy ];
sum_i += mBY[ 2 + 5*idy ];
sum_i += mBY[ 3 + 5*idy ];
mBY[ 4 + 5*idy ] = (mBY[ 4 + 5*idy ] - sum_i ) / mA_L[ 24 ] ;
}
// Solve L System DONE!
barrier(CLK_LOCAL_MEM_FENCE);
// Solve U System
// Solve U System - Row 4
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 5 + 5*idy ];
sum_i += mBY[ 4 + 5*idy ];
sum_i += mBY[ 3 + 5*idy ];
sum_i += mBY[ 2 + 5*idy ];
sum_i += mBY[ 1 + 5*idy ];
mBY[ 4 + 5*idy ] = (mBY[ 4 + 5*idy ] - sum_i ) / mA[ 24 ] ;
}
// Solve U System - Row 3
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 5 + 5*idy ];
sum_i += mBY[ 4 + 5*idy ];
sum_i += mBY[ 3 + 5*idy ];
sum_i += mBY[ 2 + 5*idy ];
mBY[ 3 + 5*idy ] = (mBY[ 3 + 5*idy ] - sum_i ) / mA[ 18 ] ;
}
// Solve U System - Row 2
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 5 + 5*idy ];
sum_i += mBY[ 4 + 5*idy ];
sum_i += mBY[ 3 + 5*idy ];
mBY[ 2 + 5*idy ] = (mBY[ 2 + 5*idy ] - sum_i ) / mA[ 12 ] ;
}
// Solve U System - Row 1
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 5 + 5*idy ];
sum_i += mBY[ 4 + 5*idy ];
mBY[ 1 + 5*idy ] = (mBY[ 1 + 5*idy ] - sum_i ) / mA[ 6 ] ;
}
// Solve U System - Row 0
if(idx==0 && idy<2 ){
float sum_i = 0.0f;
sum_i += mBY[ 5 + 5*idy ];
mBY[ 0 + 5*idy ] = (mBY[ 0 + 5*idy ] - sum_i ) / mA[ 0 ] ;
}
// Solve U System DONE!
barrier(CLK_LOCAL_MEM_FENCE);
//-------------------------------------------------------------------------------------
if(idx==0&&idy==0){
mr_sB_2_M_0 = mBY[ 6 ];
mr_sB_3_M_0 = mBY[ 7 ];
mr_sC_1_M_0 = mBY[ 8 ];
mr_sD_2_M_0 = mBY[ 9 ];
mr_sD_3_M_0 = mBY[ 10 ];
mr_sB_2_M_1 = mBY[ 11 ];
mr_sB_3_M_1 = mBY[ 12 ];
mr_sC_1_M_1 = mBY[ 13 ];
mr_sD_2_M_1 = mBY[ 14 ];
mr_sD_3_M_1 = mBY[ 15 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//-------------------------------------------------------------------------------------
//----Init A---------------------------------------------------------------------------------
if(idx==0 && idy==0){
mAm = 2;
mAn = 2;
mA[0] = ( (  -1.0f  * v_R01   )  +(  -1.0f  * v_R03   )   ) ;
mA[1] =  v_R02 ;
mA[2] =  v_R03 ;
mA[3] = ( (  -1.0f  * v_R02   )  +(  -1.0f  * v_R05   )   ) ;
}
//----Init B---------------------------------------------------------------------------
if(idx==1 && idy==0){
mBm = 2;
mBn = 2;
mB[0] = (  -1.0f  * v_R01   ) ;
mB[1] = 0.0f;
mB[2] = 0.0f;
mB[3] = (  -1.0f  * v_R05   ) ;
}
//----Init Y---------------------------------------------------------------------------
if(idx==2 && idy==0){
mYm = 2;
mYn = 3;
mY[0] = ( (  mr_sA_2_3_M_0   )   ) ;
mY[1] = ( (  mr_sB_3_M_0  * mr_sC_1_M_0   )   ) ;
mY[2] = ( (  mr_sA_2_3_M_1   )   ) ;
mY[3] = ( (  mr_sB_3_M_0  * mr_sC_1_M_1   )  +(  mr_sB_3_M_1  * mr_sC_1_M_0   )   ) ;
mY[4] = ( (  mr_sA_2_3_M_2   )   ) ;
mY[5] = ( (  mr_sB_3_M_1  * mr_sC_1_M_1   )   ) ;
}

//------COMPUTE!-----------------------------------------------------------------------
//LU-DECOMPOSITION (in place): mA_L * mA = mA
if( idx<2 && idy<2 ){ mA_L[ idx+2*idy ] = (idx==idy)?1.0f:0.0f; }
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 0 / 1
if( idx==0 && idy > 0 && idy < 2 ){
mA_L[ idy + 0 ] = mA[ idy + 0 ] / mA[ 0 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>0 && idx < 2 && idy < 2 ){
mA[ idx + 2*idy ] -= mA_L[ idx + 0 ] * mA[0+ idy*2 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 1 / 1
if( idx==0 && idy > 1 && idy < 2 ){
mA_L[ idy + 2 ] = mA[ idy + 2 ] / mA[ 3 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>1 && idx < 2 && idy < 2 ){
mA[ idx + 2*idy ] -= mA_L[ idx + 2 ] * mA[1+ idy*2 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION DONE!
//----------------------------------
//MATRIX MUL: mBY = mB * mY
//Work item (idx,idy) compute mBY[idx,idy]
if( idx < 2 && idy < 3 ) {
float mulsum = 0.0f;
mulsum += mB[ 2*0 + idx ] * mY[ idy * 2 + 0 ];
mulsum += mB[ 2*1 + idx ] * mY[ idy * 2 + 1 ];
mBY[idx+2*idy] = mulsum;
}
//MATRIX MUL COMPLETE!
barrier(CLK_LOCAL_MEM_FENCE);
// Solve L System
// Solve L System - Row 0
if(idx==0 && idy<3 ){
float sum_i = 0.0f;
mBY[ 0 + 2*idy ] = (mBY[ 0 + 2*idy ] - sum_i ) / mA_L[ 0 ] ;
}
// Solve L System - Row 1
if(idx==0 && idy<3 ){
float sum_i = 0.0f;
sum_i += mBY[ 0 + 2*idy ];
mBY[ 1 + 2*idy ] = (mBY[ 1 + 2*idy ] - sum_i ) / mA_L[ 3 ] ;
}
// Solve L System DONE!
barrier(CLK_LOCAL_MEM_FENCE);
// Solve U System
// Solve U System - Row 1
if(idx==0 && idy<3 ){
float sum_i = 0.0f;
sum_i += mBY[ 2 + 2*idy ];
sum_i += mBY[ 1 + 2*idy ];
mBY[ 1 + 2*idy ] = (mBY[ 1 + 2*idy ] - sum_i ) / mA[ 3 ] ;
}
// Solve U System - Row 0
if(idx==0 && idy<3 ){
float sum_i = 0.0f;
sum_i += mBY[ 2 + 2*idy ];
mBY[ 0 + 2*idy ] = (mBY[ 0 + 2*idy ] - sum_i ) / mA[ 0 ] ;
}
// Solve U System DONE!
barrier(CLK_LOCAL_MEM_FENCE);
//-------------------------------------------------------------------------------------
if(idx==0&&idy==0){
mr_sB_2_3_M_0 = mBY[ 3 ];
mr_sD_2_3_M_0 = mBY[ 4 ];
mr_sB_2_3_M_1 = mBY[ 5 ];
mr_sD_2_3_M_1 = mBY[ 6 ];
mr_sB_2_3_M_2 = mBY[ 7 ];
mr_sD_2_3_M_2 = mBY[ 8 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//-------------------------------------------------------------------------------------
//----Init A---------------------------------------------------------------------------------
if(idx==0 && idy==0){
mAm = 3;
mAn = 3;
mA[0] = ( (  -1.0f  * v_R03   )  +(  -1.0f  * v_R01   )   ) ;
mA[1] =  v_R02 ;
mA[2] = 0.0f;
mA[3] =  v_R03 ;
mA[4] = ( (  -1.0f  * v_R05   )  +(  -1.0f  * v_R02   )   ) ;
mA[5] =  v_R06 ;
mA[6] = 0.0f;
mA[7] = 0.0f;
mA[8] = ( (  -1.0f  * v_R06   )   ) ;
}
//----Init B---------------------------------------------------------------------------
if(idx==1 && idy==0){
mBm = 3;
mBn = 2;
mB[0] = 0.0f;
mB[1] = (  -1.0f  * v_R05   ) ;
mB[2] = 0.0f;
mB[3] = (  -1.0f  * v_R01   ) ;
mB[4] = 0.0f;
mB[5] = 0.0f;
}
//----Init Y---------------------------------------------------------------------------
if(idx==2 && idy==0){
mYm = 2;
mYn = 4;
mY[0] = ( (  mr_sB_2_3_M_0  * mr_sC_1_M_0   )   ) ;
mY[1] = ( (  mr_sA_1_2_3_M_0   )   ) ;
mY[2] = ( (  mr_sB_2_3_M_0  * mr_sC_1_M_1   )  +(  mr_sB_2_3_M_1  * mr_sC_1_M_0   )   ) ;
mY[3] = ( (  mr_sA_1_2_3_M_1   )   ) ;
mY[4] = ( (  mr_sB_2_3_M_1  * mr_sC_1_M_1   )  +(  mr_sB_2_3_M_2  * mr_sC_1_M_0   )   ) ;
mY[5] = ( (  mr_sA_1_2_3_M_2   )   ) ;
mY[6] = ( (  mr_sB_2_3_M_2  * mr_sC_1_M_1   )   ) ;
mY[7] = ( (  mr_sA_1_2_3_M_3   )   ) ;
}

//------COMPUTE!-----------------------------------------------------------------------
//LU-DECOMPOSITION (in place): mA_L * mA = mA
if( idx<3 && idy<3 ){ mA_L[ idx+3*idy ] = (idx==idy)?1.0f:0.0f; }
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 0 / 2
if( idx==0 && idy > 0 && idy < 3 ){
mA_L[ idy + 0 ] = mA[ idy + 0 ] / mA[ 0 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>0 && idx < 3 && idy < 3 ){
mA[ idx + 3*idy ] -= mA_L[ idx + 0 ] * mA[0+ idy*3 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 1 / 2
if( idx==0 && idy > 1 && idy < 3 ){
mA_L[ idy + 3 ] = mA[ idy + 3 ] / mA[ 4 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>1 && idx < 3 && idy < 3 ){
mA[ idx + 3*idy ] -= mA_L[ idx + 3 ] * mA[1+ idy*3 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION STEP 2 / 2
if( idx==0 && idy > 2 && idy < 3 ){
mA_L[ idy + 6 ] = mA[ idy + 6 ] / mA[ 8 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
if( idx>2 && idx < 3 && idy < 3 ){
mA[ idx + 3*idy ] -= mA_L[ idx + 6 ] * mA[2+ idy*3 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//LU-DECOMPOSITION DONE!
//----------------------------------
//MATRIX MUL: mBY = mB * mY
//Work item (idx,idy) compute mBY[idx,idy]
if( idx < 3 && idy < 4 ) {
float mulsum = 0.0f;
mulsum += mB[ 3*0 + idx ] * mY[ idy * 2 + 0 ];
mulsum += mB[ 3*1 + idx ] * mY[ idy * 2 + 1 ];
mBY[idx+3*idy] = mulsum;
}
//MATRIX MUL COMPLETE!
barrier(CLK_LOCAL_MEM_FENCE);
// Solve L System
// Solve L System - Row 0
if(idx==0 && idy<4 ){
float sum_i = 0.0f;
mBY[ 0 + 3*idy ] = (mBY[ 0 + 3*idy ] - sum_i ) / mA_L[ 0 ] ;
}
// Solve L System - Row 1
if(idx==0 && idy<4 ){
float sum_i = 0.0f;
sum_i += mBY[ 0 + 3*idy ];
mBY[ 1 + 3*idy ] = (mBY[ 1 + 3*idy ] - sum_i ) / mA_L[ 4 ] ;
}
// Solve L System - Row 2
if(idx==0 && idy<4 ){
float sum_i = 0.0f;
sum_i += mBY[ 0 + 3*idy ];
sum_i += mBY[ 1 + 3*idy ];
mBY[ 2 + 3*idy ] = (mBY[ 2 + 3*idy ] - sum_i ) / mA_L[ 8 ] ;
}
// Solve L System DONE!
barrier(CLK_LOCAL_MEM_FENCE);
// Solve U System
// Solve U System - Row 2
if(idx==0 && idy<4 ){
float sum_i = 0.0f;
sum_i += mBY[ 3 + 3*idy ];
sum_i += mBY[ 2 + 3*idy ];
sum_i += mBY[ 1 + 3*idy ];
mBY[ 2 + 3*idy ] = (mBY[ 2 + 3*idy ] - sum_i ) / mA[ 8 ] ;
}
// Solve U System - Row 1
if(idx==0 && idy<4 ){
float sum_i = 0.0f;
sum_i += mBY[ 3 + 3*idy ];
sum_i += mBY[ 2 + 3*idy ];
mBY[ 1 + 3*idy ] = (mBY[ 1 + 3*idy ] - sum_i ) / mA[ 4 ] ;
}
// Solve U System - Row 0
if(idx==0 && idy<4 ){
float sum_i = 0.0f;
sum_i += mBY[ 3 + 3*idy ];
mBY[ 0 + 3*idy ] = (mBY[ 0 + 3*idy ] - sum_i ) / mA[ 0 ] ;
}
// Solve U System DONE!
barrier(CLK_LOCAL_MEM_FENCE);
//-------------------------------------------------------------------------------------
if(idx==0&&idy==0){
mr_sB_1_2_3_M_0 = mBY[ 4 ];
mr_sD_1_2_3_M_0 = mBY[ 5 ];
mr_sF_1_2_3_M_0 = mBY[ 6 ];
mr_sB_1_2_3_M_1 = mBY[ 7 ];
mr_sD_1_2_3_M_1 = mBY[ 8 ];
mr_sF_1_2_3_M_1 = mBY[ 9 ];
mr_sB_1_2_3_M_2 = mBY[ 10 ];
mr_sD_1_2_3_M_2 = mBY[ 11 ];
mr_sF_1_2_3_M_2 = mBY[ 12 ];
mr_sB_1_2_3_M_3 = mBY[ 13 ];
mr_sD_1_2_3_M_3 = mBY[ 14 ];
mr_sF_1_2_3_M_3 = mBY[ 15 ];
}
barrier(CLK_LOCAL_MEM_FENCE);
//-------------------------------------------------------------------------------------
//---RESULTS-----------------------------------------------------------------------------------------
if(idx==0&&idy==0){
mds_out[ offset_target + 1 ] = mr_sF_1_2_3_M_0;
mds_out[ offset_target + 2 ] = mr_sF_1_2_3_M_1;
mds_out[ offset_target + 3 ] = mr_sF_1_2_3_M_2;
mds_out[ offset_target + 4 ] = mr_sF_1_2_3_M_3;
}
//---------------------------------------------------------------------------------------------------
}
"


function test_0()
    a = convert( Array{Float32,1} ,[100.; 110 ; 50 ; 20 ; 20 ; 80 ; 0;0;0;0;0;0;0;0;0;0] )
    b = convert( Array{Float32,1} ,[0.00 ; 1.00; 1.00; 0.00; 0.00; 1.0 ;0.0;0.0;1.0;0.0;0.0;0;0;0;0;0] )

    device, ctx, queue = OpenCL.cl.create_compute_context()

    a_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=a)
    b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=b)
    c_buff = cl.Buffer(Float32, ctx, :w, length(a))

    kernel_str = str_kernel_test_01
    p = cl.Program(ctx, source=kernel_str) |> cl.build!

    #p = cl.Program(ctx, source=kernel_02) |> cl.build!

    #cl.CL_device_info(1)

    k = cl.Kernel(p, "fwsim")
    queue(k, [16;16], [16;16], a_buff, b_buff, c_buff)

    #cl.api.clGetKernelWorkGroupInfo(cl.CL_KERNEL_WORK_GROUP_SIZE)
    r = cl.read(queue, c_buff)

    cl.work_group_info( k , cl.CL_KERNEL_WORK_GROUP_SIZE ,device )
    cl.work_group_info( k , cl.CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE  ,device )
end

test_0()
