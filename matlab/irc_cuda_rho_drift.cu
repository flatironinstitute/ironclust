/**
 * jrc3_cuda_rho_drift.cu
 * J. James Jun, Flatiron Institute, 2018 Jun 20
 * All sites simultaneously
*/

#include <cuda_runtime.h>
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define NC 45 //max dimm
#define CHUNK 16
#define SINGLE_INF (3.402E+38) // equipvalent to NAN. consider -1 value


__host__ __device__ float calc_dist_(const float mrFet_[][NC], const int i_c, const float vrFet_tx[], const int nC){
    float dist = 0.0f;    
    for (int iC=0; iC<nC; ++iC){
        float temp = mrFet_[i_c][iC] - vrFet_tx[iC];
        dist += temp * temp;
    }    
    return dist;
}


__global__ void jrc3_cuda_rho_drift(float * vrRho1, const int * miSite12, const float * mrFet1, const float * mrFet2, const int * vnConst, const float dc2){
    int i1 = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK;   // base index of i1    
    int tx = threadIdx.x; //nThreads for i12 index    
    int i1_tx = i1 + tx;
    int n1 = vnConst[0];
    int n12 = vnConst[1];
    int nC = vnConst[2];
    
    __shared__ float mrFet1_[CHUNK][NC]; //, mrFet2_[CHUNK][NC];
    __shared__ int viSite1_[CHUNK]; //, viSite2_[CHUNK];    
    
    __shared__ int mnRho1_[NTHREADS][CHUNK]; // count then divide later    
    __shared__ int mnComp1_[NTHREADS][CHUNK]; // count number of elements compared
    
    // cache shared memory
    if (tx < nC){ //use tx as iC
        for (int i_c = 0; i_c < CHUNK; ++i_c){
            int i1_c = i_c + i1;
            if (i1_c < n1){
                mrFet1_[i_c][tx] = mrFet1[tx + i1_c * nC];
                //mrFet2_[i_c][tx] = mrFet2[tx + i1_c * nC];
                viSite1_[i_c] = miSite12[i1_c + 0 * n12];
                //viSite2_[i_c] = miSite12[i1_c + 1 * n12];
            }else{
                mrFet1_[i_c][tx] = 0.0f;
                //mrFet2_[i_c][tx] = 0.0f;
                viSite1_[i_c] = -1;
                //viSite2_[i_c] = -1;
            }
        }        
    }    
    for (int i_c = 0; i_c < CHUNK; ++i_c){        
        mnRho1_[tx][i_c] = 0; // initialize rho
        mnComp1_[tx][i_c] = 0;        
    }     
    __syncthreads();        

    
    // Inspect distance relationship between i1 and i12_tx
    for (int i12_tx = tx; i12_tx < n12; i12_tx += blockDim.x){        
        // Cache info
        int iSite1_tx = miSite12[i12_tx + 0 * n12]; // primary site
        int iSite2_tx = miSite12[i12_tx + 1 * n12]; // secondary site
        float vrFet1_tx[NC], vrFet2_tx[NC];
        //char fFet1 = 0, fFet2 = 0; // cache status
        for (int iC = 0; iC < nC; ++iC){
            vrFet1_tx[iC] = mrFet1[iC + i12_tx * nC];
            vrFet2_tx[iC] = mrFet2[iC + i12_tx * nC];           
        }
        
        //Compute distance and rho
        for (int i_c = 0; i_c < CHUNK; ++i_c){
            float dist = 0.0f;
            if (viSite1_[i_c] == iSite1_tx){
                dist = calc_dist_(mrFet1_, i_c, vrFet1_tx, nC);
            }else if (viSite1_[i_c] == iSite2_tx){
                dist = calc_dist_(mrFet1_, i_c, vrFet2_tx, nC);
            }
            else{
                continue;
            }
            
            if (dist <= dc2) ++mnRho1_[tx][i_c];
            ++mnComp1_[tx][i_c];
        }        
    } // while
    
    // final count
    __syncthreads();
    if (tx < CHUNK){  // use tx as i_c
        int nRho1 = 0;
        int nComp1 = 0;
        for (int tx1=0; tx1<blockDim.x; ++tx1){
            nRho1 += mnRho1_[tx1][tx];
            nComp1 += mnComp1_[tx1][tx];
        }
        if (i1_tx < n1){
            vrRho1[i1_tx] = (float)(((double)(nRho1)) / ((double)nComp1));
        }
    }
} // func