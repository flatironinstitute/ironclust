/* 
 * James Jun 2018/7/3 
 * Fast approximation of knn using binned minimum parallel search
 * J. James Jun, Flatiron Institute, 2018 July 5
*/

#include <cuda_runtime.h>
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? -(my_val) : (my_val)
#define NC (45) //3pca x 16 channels max
#define SINGLE_INF (3.402E+38) // equipvalent to NAN. consider -1 value
#define SWAP(x, y, T) do { T SWAP = x; x = y; y = SWAP; } while (0)
#define CHUNK (4) // 4
#define NTHREADS (512) // 512
#define NM (256) // number of drift
//#define CHUNK (8)
//#define NTHREADS (256)

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
__global__ void cuda_knn_drift(float *K, int *I, float const *B, float const *A, unsigned char const *MB, unsigned char const *MA, char const *M, const int *vnConst){
    
    int nB = vnConst[0];
    int nA = vnConst[1];
    int nC = vnConst[2];
    int knn = vnConst[3]; // not used for now
    int nM = vnConst[4]; // nD <= ND

    int tx = threadIdx.x;
    int iA = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK + tx;
    int nThreads = blockDim.x; // must be less than NTHREADS

    __shared__ float sA[NC][CHUNK];
    __shared__ float sD[NTHREADS][CHUNK];
    __shared__ int sI[NTHREADS][CHUNK];
    __shared__ char sM[NM][CHUNK]; // logical
    
    // initialize
    if (tx < CHUNK){        
        if (iA < nA){
            int iA_ = tx; // loop over CHUNK
            for (int iC=0; iC<nC; ++iC) sA[iC][iA_] = A[iC + iA*nC]; // copy A->sA   
            for (int iT=0; iT<nThreads; ++iT){
                sD[iT][iA_] = SINGLE_INF; // sD = inf
                sI[iT][iA_] = 0;                
            }
            for (int iM=0; iM<nM; ++iM){
                int ma_ = (int)MA[iA];
                sM[iM][iA_] = M[iM + ma_ * nM];
            }
        }        
    }
    __syncthreads();
    
    
    // find minimum for each bin, stride of 2xnThread to save shared memory
    for (int iB=tx; iB<nB; iB+=nThreads){    
        // Initialize distance vector
        float dist_[CHUNK];     
        for (int iA_=0; iA_<CHUNK; ++iA_) dist_[iA_] = 0.0f;
        
        int mb_ = (int)MB[iB];
        int n_inf_ = 0;
        for (int iA_=0; iA_<CHUNK; ++iA_){            
            if (sM[mb_][iA_]==0){
                dist_[iA_] = SINGLE_INF;
                n_inf_++;
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
        int iB_ = tx;
        for (int iA_=0; iA_<CHUNK; ++iA_){
            float d_ = dist_[iA_];
            if (dist_[iA_] < sD[iB_][iA_] && d_>0){
                sD[iB_][iA_] = d_;
                sI[iB_][iA_] = iB;
            }
        }
    } // while
    __syncthreads();
    
    
    // sort up to kth element using bubble sort
    if (tx < CHUNK){           
        if (iA < nA){
            int iA_ = tx;
            for (int iK=0; iK<knn; ++iK){
                float dmin_ = sD[iK][iA_];
                int imin_ = iK;
                for (int iB_=iK; iB_<nThreads; ++iB_){
                    float d_ = sD[iB_][iA_];
                    if (d_ < dmin_){
                        dmin_ = d_;
                        imin_ = iB_;
                    }
                }
                SWAP(sD[imin_][iA_], sD[iK][iA_], float);
                SWAP(sI[imin_][iA_], sI[iK][iA_], int);
            }      
            // average or take the last?
            K[iA] = sqrt(ABS(sD[knn-1][iA_]));
            for (int iK=0; iK<knn; ++iK){
                I[iK + knn*iA] = sI[iK][iA_] + 1; // matlab 1 base
            }
        }        
    } // if
} // func
