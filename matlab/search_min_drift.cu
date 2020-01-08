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
#define CHUNK (8) // 16

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * D: Mininum distances (nT x nA, float)
 * I: Index (nT x nA, int32)
 * F: Feature matrix (nC x nF, single)
 * IBB: vector of B offset (nBB x 1, int)
 * NBB: vector of B counts (nBB x 1, int)
 * const: [nC, nF, nBB, iA0, nA]
 */
__global__ void search_min_drift(float *D, int *I, float const *F, const int *IBB, const int *NBB, const int *vnConst){
    
    int nC = vnConst[0]; // number of feature dimension
    int nF = vnConst[1]; // number of channels
    int nBB = vnConst[2]; // number of blocks to read
    int iA0 = vnConst[3]; // A offset
    int nA = vnConst[4]; // number of A to read

    int tx = threadIdx.x;
    int nT = blockDim.x; // must be less than NTHREADS

    // shared memory
    __shared__ float sA[NC][CHUNK];
    __shared__ int sI[CHUNK];
    // thread memory
    float tD[CHUNK]; // t: thread, s: shraed
    int tI[CHUNK];
    
    // initialize
    if (tx < CHUNK){    
        int iA_ = tx; // loop over CHUNK
        int iA = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK + tx;
        iA = (iA % nA) + iA0;
        for (int iC=0; iC<nC; ++iC){
            sA[iC][iA_] = F[iC + iA*nC]; // copy A->sA 
        }
        sI[iA_] = iA;         
    }    
    __syncthreads();
    
    //#pragma unroll
    for (int iA_=0; iA_<CHUNK; ++iA_){            
        tD[iA_] = SINGLE_INF; // sD = inf
        tI[iA_] = sI[iA_];                
    }     
    
    // find minimum for each bin, stride of 2xnThread to save shared memory
    int iB0 = IBB[tx % nBB];
    int nB0 = NBB[tx % nBB];
    int iT0 = tx / nBB;
    int nT0 = nT / nBB;
    for (int iB=iT0+iB0; iB<nB0+iB0; iB+=nT0){    
        // Initialize distance vector
        float dist_[CHUNK];  // #programa unroll?   
        for (int iA_=0; iA_<CHUNK; ++iA_) dist_[iA_] = 0.0f;

        for (int iC=0; iC<nC; ++iC){
            float b_ = F[iC + iB*nC];
            for (int iA_=0; iA_<CHUNK; ++iA_){
                float d_ = b_ - sA[iC][iA_];
                dist_[iA_] += (d_ * d_);
            }            
        }         
        for (int iA_=0; iA_<CHUNK; ++iA_){
            float d_ = dist_[iA_];
            if (dist_[iA_] < tD[iA_] && d_ > 0.0f){
                tD[iA_] = d_;
                tI[iA_] = iB;
            }
        }
    } // while
    
    // write the output
    for (int iA_=0; iA_<CHUNK; ++iA_){ 
        int iA1 = sI[iA_] - iA0; // guaranteed to be bound due to modulus algebra     
        if (iA1 >= 0 && iA1 < nA){
            D[tx + nT*iA1] = sqrt(ABS(tD[iA_]));
            I[tx + nT*iA1] = tI[iA_] + 1; // matlab 1 base
        }
    } // for
} // func
