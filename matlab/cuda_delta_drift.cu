/**
 * cuda_delta_knn.cu
 * block loading delta calculation.
 * system('nvcc -ptx -m 64 -arch sm_35 jrc3_cuda_rho.cu')
 * J. James Jun, Flatiron Institute, 2018 July 5
*/

#include <cuda_runtime.h>
// #include "cublas_v2.h"
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS (128)
#define NC (45) // number of Channels
#define CHUNK (16) //previously defined as CHUNK
#define SINGLE_INF (3.402E+38)
#define NM (256) // number of drift

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * K: K nearest neighbor
 * I: Index (knn x nA, int32)
 * B: Feature matrix (nC x nB, single)
 * A: Feature matrix (nC x nB, single)
 * MB: Drift index B (nB x 1, unit8)
 * MA: Drift Index A (nA x 1, uint8)
 * M: drift connection matrix (nD x nD, logical)
 */

__global__ void cuda_delta_drift(float *D, unsigned int *N, const float *B, const float *A, unsigned char const *MB, unsigned char const *MA, char const *M, const float *R_B, const float *R_A, const int *vnConst){
    
    int nB = vnConst[0];
    int nA = vnConst[1];
    int nC = vnConst[2];
    int nM = vnConst[3]; // nD <= ND
    
    int tx = threadIdx.x;
    int iA = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK + tx;
    int nThreads = blockDim.x; // must be less than NTHREADS

    __shared__ float sA[NC][CHUNK], sR_A[CHUNK]; // march through A
    __shared__ float sD[NTHREADS][CHUNK]; // count then divide later    
    __shared__ unsigned int sN[NTHREADS][CHUNK]; 
    __shared__ char sM[NM][CHUNK]; // logical
    
    // initialize
    if (tx < CHUNK){         
        if (iA < nA){
            int iA_ = tx;
            for (int iC=0; iC<nC; ++iC) sA[iC][tx] = A[iC + iA*nC]; // copy A->sA
            sR_A[iA_] = R_A[iA]; // copy R_A->sR_A
            for (int iB_=0; iB_<nThreads; ++iB_) sD[iB_][iA_] = SINGLE_INF; // sD = inf        
            for (int iM=0; iM<nM; ++iM){
                int ma_ = (int)MA[iA];
                sM[iM][iA_] = M[iM + ma_ * nM];
            }        
        }
    }
    __syncthreads();        

    
    // Search min distance having a greater rho
    for (int iB=tx; iB<nB; iB+=nThreads){        
        // Initialize distance vector
        float dist_[CHUNK];        
        for (int iA_=0; iA_<CHUNK; ++iA_) dist_[iA_] = 0.0f;        
        
        int mb_ = (int)MB[iB];
        float rb_ = R_B[iB];
        int n_inf_ = 0;
        for (int iA_=0; iA_<CHUNK; ++iA_){            
            if (sM[mb_][iA_]==0 || rb_ <= sR_A[iA_]){
                dist_[iA_] = SINGLE_INF;
                n_inf_++;
            }
        }         
        if (n_inf_ >= CHUNK) continue; 
        
        for (int iC=0; iC<nC; ++iC){
            float b_ = B[iC + iB * nC];
            for (int iA_=0; iA_<CHUNK; ++iA_){
                float d_ = b_ - sA[iC][iA_];
                dist_[iA_] += (d_ * d_);
            }            
        }
        int iB_ = tx;
        for (int iA_=0; iA_<CHUNK; ++iA_){   
            if (dist_[iA_] < sD[iB_][iA_]){
                sD[iB_][iA_] = dist_[iA_];
                sN[iB_][iA_] = iB;
            }
        }
    } // while    
    __syncthreads();    
    
    
    // final count
    if (tx < CHUNK){        
        if (iA < nA){
            int iA_ = tx;
            float dmin_ = SINGLE_INF;
            unsigned int imin_ = iA; // point to self initially
            for (int iB_=0; iB_<nThreads; ++iB_){                
                if (sD[iB_][iA_] < dmin_){
                    dmin_ = sD[iB_][iA_];
                    imin_ = sN[iB_][iA_];
                }
            }
            D[iA] = sqrtf(ABS(dmin_));
            N[iA] = imin_ + 1; //Matlab index output
        }
    }
} // func