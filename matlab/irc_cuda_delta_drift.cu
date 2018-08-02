/**
 * jrc3_cuda_delta_drift.cu
 * J. James Jun, Flatiron Institute, 2018 Jun 20
 * All sites simultaneously
*/

#include <cuda_runtime.h>
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define NC 45 // number of Channels
#define CHUNK 16 //previously defined as CHUNK
#define SINGLE_INF (3.402E+38)


__host__ __device__ float calc_dist_(const float mrFet_[][NC], const int i_c, const float vrFet_tx[], const int nC){
    float dist = 0.0f;    
    for (int iC=0; iC<nC; ++iC){
        float temp = mrFet_[i_c][iC] - vrFet_tx[iC];
        dist += temp * temp;
    }    
    return dist;
}


__global__ void jrc3_cuda_delta_drift(float * vrDelta1, unsigned int * viNneigh1, const int * miSite12, const float * mrFet1, const float * mrFet2, const float * vrRho1, const int * vnConst, const float dc2){
    int i1 = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK;   // base index of i1    
    int tx = threadIdx.x; //nThreads for i12 index    
    int i1_tx = i1 + tx;
    int n1 = vnConst[0];
    int n12 = vnConst[1];
    int nC = vnConst[2];
    
    __shared__ float mrFet1_[CHUNK][NC]; //, mrFet2_[CHUNK][NC];
    __shared__ int viSite1_[CHUNK]; //, viSite2_[CHUNK];    
    
    __shared__ float vrRho1_[CHUNK];
    __shared__ float mrDelta1_[NTHREADS][CHUNK];
    __shared__ unsigned int miNneigh1_[NTHREADS][CHUNK]; 
    
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
    if (tx < CHUNK && i1_tx < n1){
        vrRho1_[tx] = vrRho1[i1_tx];
    }

    float vr_minDist1[CHUNK];
    unsigned int vi_minIdx1[CHUNK];
    for (int i_c = 0; i_c < CHUNK; ++i_c){
        vr_minDist1[i_c] = SINGLE_INF;
        vi_minIdx1[i_c] = i1 + i_c; // self
    }        
    __syncthreads();  
    
    
    // fill in the shared memory A
    for (int i12_tx = tx; i12_tx < n12; i12_tx += blockDim.x){
        float rho1_tx = vrRho1[i12_tx];

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
            if (rho1_tx <= vrRho1_[i_c]) continue;
            
            float dist = 0.0f;
            if (viSite1_[i_c] == iSite1_tx){
                dist = calc_dist_(mrFet1_, i_c, vrFet1_tx, nC);
            }else if (viSite1_[i_c] == iSite2_tx){
                dist = calc_dist_(mrFet1_, i_c, vrFet2_tx, nC);
            }
            else{
                continue;
            }
            
            // compare if its the min dist
            if (dist < vr_minDist1[i_c]){
                vr_minDist1[i_c] = dist;
                vi_minIdx1[i_c] = i12_tx;
            }
        }
    } // while
    
    // collect result from each thread
    for (int i_c = 0; i_c < CHUNK; ++i_c){        
        mrDelta1_[tx][i_c] = vr_minDist1[i_c];
        miNneigh1_[tx][i_c] = vi_minIdx1[i_c];
    }
    __syncthreads();    
    
    // final count    
    if (tx < CHUNK){
        float minDist1 = SINGLE_INF;
        unsigned int minIdx1 = i1_tx;
        for (int tx1=0; tx1<blockDim.x; ++tx1){
            if (mrDelta1_[tx1][tx] < minDist1){
                minDist1 = mrDelta1_[tx1][tx];
                minIdx1 = miNneigh1_[tx1][tx];
            }
        }
        if (i1_tx < n1){
            vrDelta1[i1_tx] = sqrtf(ABS(minDist1) / dc2); 
            viNneigh1[i1_tx] = minIdx1 + 1; //Matlab index output
        }
    }
} // func