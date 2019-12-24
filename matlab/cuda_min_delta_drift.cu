/* 
 * James Jun 2019/12/22 
 * Fast approximation of knn using binned minimum parallel search
*/

#include <cuda_runtime.h>
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? -(my_val) : (my_val)
#define NC (45) //3pca x 16 channels max
#define SINGLE_INF (3.402E+38) // equipvalent to NAN. consider -1 value
#define SWAP(x, y, T) do { T SWAP = x; x = y; y = SWAP; } while (0)
#define CHUNK (32) // 4
#define NM (256) // number of drift

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * D: Mininum distances (nT x nA, float)
 * I: Index (nT x nA, int32)
 * B: Feature matrix (nC x nB, single)
 * A: Feature matrix (nC x nB, single)
 * MB: Drift index B (nB x 1, unit8)
 * MA: Drift Index A (nA x 1, uint8)
 * M: drift connection matrix (nD x nD, logical)
 * const: [nB, nA, nC, nT, nM, int32], nT: # threads
 */
__global__ void cuda_min_delta_drift(float *D, int *I, float const *B, float const *A, unsigned char const *MB, unsigned char const *MA, char const *M, const float *R_B, const float *R_A, const int *vnConst){
    
    int nB = vnConst[0];
    int nA = vnConst[1];
    int nC = vnConst[2];
    int nM = vnConst[3]; // nD <= ND

    int tx = threadIdx.x;
    int nT = blockDim.x; // must be less than NTHREADS

    // shared memory
    __shared__ float sA[NC][CHUNK], sR_A[CHUNK];
    __shared__ int sI[CHUNK];
    __shared__ char sM[NM][CHUNK]; // logical
    // thread memory
    float tD[CHUNK]; // t: thread, s: shraed
    int tI[CHUNK];
    
    // initialize
    if (tx < CHUNK){    
        int iA_ = tx; // loop over CHUNK
        int iA = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK + tx;
        iA = iA % nA;
        sR_A[iA_] = R_A[iA]; // copy R_A->sR_A
        for (int iC=0; iC<nC; ++iC){
            sA[iC][iA_] = A[iC + iA*nC]; // copy A->sA 
        }
        sI[iA_] = iA;
        for (int iM=0; iM<nM; ++iM){
            int ma_ = (int)MA[iA];
            sM[iM][iA_] = M[iM + ma_ * nM];
        }           
    }    
    __syncthreads();
    
    //#pragma unroll
    for (int iA_=0; iA_<CHUNK; ++iA_){            
        tD[iA_] = SINGLE_INF; // sD = inf
        tI[iA_] = sI[iA_];                
    }     
    
    // find minimum for each bin, stride of 2xnThread to save shared memory
    for (int iB=tx; iB<nB; iB+=nT){    
        // Initialize distance vector
        float dist_[CHUNK];  // #programa unroll?   
        
        int mb_ = (int)MB[iB];
        float rb_ = R_B[iB];
        int n_inf_ = 0;
        for (int iA_=0; iA_<CHUNK; ++iA_){            
            if (sM[mb_][iA_]==0 || rb_ <= sR_A[iA_]){
                dist_[iA_] = SINGLE_INF;
                n_inf_++;
            }else{
                dist_[iA_] = 0.0f;
            }
        }        
        if (n_inf_ >= CHUNK) continue;

        for (int iC=0; iC<nC; ++iC){
            float b_ = B[iC + iB*nC];
            for (int iA_=0; iA_<CHUNK; ++iA_){
                float d_ = b_ - sA[iC][iA_];
                dist_[iA_] += (d_ * d_);
            }            
        }         
        for (int iA_=0; iA_<CHUNK; ++iA_){
            float d_ = dist_[iA_];
            if (dist_[iA_] < tD[iA_]){
                tD[iA_] = d_;
                tI[iA_] = iB;
            }
        }
    } // while
    
    // write the output
    for (int iA_=0; iA_<CHUNK; ++iA_){ 
        int iA = sI[iA_];
        if (iA<nA && iA >= 0){
            D[tx + nT*iA] = sqrt(ABS(tD[iA_]));
            I[tx + nT*iA] = tI[iA_] + 1; // matlab 1 base
        }
    } // for
} // func
