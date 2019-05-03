extern "C"{
#include "cluster.h"
}
#include "utils.h"
#include <cuda.h>
#include <cuda_runtime.h>

__global__ void update_colcand(bool *dcolcand, discrete *dg1, discrete *dg2, int dcols)
{
	int tx = threadIdx.x;
	for(int i = 0; i < tx; ++i)
	{	
		if((i < dcols))
			if(dcolcand[i] && (dg1[i] != dg2[i]))
				dcolcand[i] = FALSE;
	}
}

__global__ void profile_init(int *dprofiles, int *drows)
{
	int tx = threadIdx.x;
	for(int i = 0; i < tx; ++i)
		if(i < drows)
			dprofiles[i] = 0;
}

__global__ void clusprofile_init(int dcols, int dsigman, int *dprofile)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	for(int i = 0; i < tx; ++i)
	{
		if(i < dcols)
		{
			for(int k = 0; k < ty; ++i)
			{
				if(k < dsigma)
					dprofile[i][k] = 0;
			}
		}
	}
}

__global__ void colcand_init(bool *dcolcand, int dcols)
{
	int tx = threadIdx.x;
	for(int i = 0; i < tx; ++i)
		if(i < dcols)
			dcolcand[i] = FALSE
}

__global__ void candidates_init(bool *dcandidates, int drows)
{
	int tx = threadIdx.x;
	for(int i = 0; i < tx; ++i)
		if(i < drows)
			candidates[i] = TRUE;
}
__global__ void intersect_row(bool *dcolcand, discrete *dg1, discrete *dg2, int dcols, int dcnt)
{
	int tx = threadIdx.x;
	int cnt = 0;
	for(int i = 0; i < tx; ++i)
	{
		if((i < dcols) && dcolcand[i] && (dg1[i] == dg2[i]) && (dg1[i] != 0))
			cnt++;
	}
	__syncthreads();
	dcnt = cnt;
}

__global__ static void reverse_row(const bool *dcolcand, discrete *dg1, discrete *dg2, int dcols, int dcnt, discrete *dsymbols)
{
	int tx = threadIdx.x;
	int cnt = 0;
	for(int i = 0; i < tx; ++i)
	{
		if((i < dcols) && dcolcand [i] && (dsymbols[dg1[i]] == -dsymbols[dg2[i]]))
			cnt++;
	}
	__syncthreads();
	dcnt = cnt;
}

__global__ static void seed_current_modify(const discrete *ds, bool *dcolcand, int *dcnt, 
						int dcomponents, int dcols, int dsigma, double dtolerance, int **dprofile)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	int indy = blockDim.y * blockIdx.y + threadIdx.y;
	int k, flag, n;
	int threshold = (dcomponents * dtolerance);
	discrete ss;
	*dcnt = 0;
	for(int i = 0; i < indx; ++i)
	{
		if( i < dcols)
		{
			flag = 0; ss = ds[i];
			for( k = 1; k < indy; ++k)
			{
				if(k < dsigma)
				{
					n = dprofile[i][k];
					if(k == ss)
						++n;
					if( n = threshold)
					{
						flag = k;
						break;
					}
				}
			}
			__syncthreads();
			if(flag)
			{
				*(dcnt)++;
				dcolcand[i] = TRUE;
			}
		}
	}
}

__global__ void bi_colcand_init(bool *dcolcand, int dcols, discrete *dg1, discrete *dg2)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	int indy = blockDim.y * blockIdx.y + threadIdx.y;
	int i;
	for(i = 0; i < indx; ++i)
		if ( i < dcols)
			dcolcand[i] = FALSE;
	for(i = 0; i < indy; ++i)
		if(i < dcols)
			dcolcand[i] = TRUE;
}

__global__ void bi_candidate_100row(bool *dcandidates, int dtop, int *darr_rows)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	for(int i = 0; i < indx; ++i)
		if(darr_rows[i] < dtop)
			dcandidates[i] = FALSE;
}

__global__ void profile_increment(Block **dbb, int *db, int *dprofiles)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	for(int i = 0; i < indx; ++i)
		if(i < dbb[db]->block_rows)
			dprofiles[dsItem(dbb[db]->genes, i)]++;
}

__global__ void get_pvalue(continuous da, int db, long double dpval)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	int i = 0;
	long double pvalue = 0;
	long double poisson = 1/exp(da);
	for(i; i < indx; ++i)
	{
		if(i > (db-1))
			pvalue = pvalue + poisson;
		else
			poisson = poisson * da/(i+1);
	}
	__syncthreads();
	dpval = pvalue;
}

void bi_colcand_wrap(bool *dcolcand, int dcols, discrete *dg1, discrete *dg2)
{
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        bi_colcand_init<<<BLOCK, GRID>>>(dcolcand, dcols, dg1, dg2);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
}
void bi_candidate_wrap(bool *dcandidates, int dtop, int *darr_rows)
{
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        bi_candidate_100row<<<BLOCK, GRID>>>(dcandidates, dtop, darr_rows);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
}
void d_profile_increment(Block **dbb, int *db1, int *db2, int *dprofiles)
{
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        profile_increment<<<BLOCK, GRID>>>(dbb, db1, dprofiles);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
        profile_increment<<<BLOCK, GRID>>>(dbb, db2, dprofiles);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
}

void d_profile_init(int *dprofiles, int *drows)
{
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        profile_init<<<BLOCK, GRID>>>(dprofiles, drows);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
}

void d_colcand_init(bool *dcolcand, int dcols)
{
	dim3 BLOCK(BLOCKsize,1,1);
	dim3 GRID(cols/BLOCKsize,1,1);
	colcand_init<<<BLOCK,GRID>>>(dcolcand,dcols);
	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());
}

void d_candidates_init(bool *dcandidates, int drows)
{
	dim3 BLOCK(BLOCKsize, 1,1);
	dim3 GRID(rows/BLOCKsize,1,1);
	candidates_init<<<BLOCK,GRID>>>(dcandidates, drows);
	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());
}

void d_seed_wrap(const discrete *ds, bool *dcolcand, int *dcnt, int dcomponents, int dcols, int dsigma, double dtolerance, int **dprofile)
{
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        seed_current_modify<<<BLOCK, GRID>>>(ds, dcolcand, dcnt, dthreshold, dcomponents, dcols, dsigma, dtolerance, dprofile);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
}

void d_reverse_wrap(const bool *dcolcand, discrete *dg1, discrete *dg2, int dcols, int dcnt, discrete *dsymbols)
{
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        reverse_row<<<BLOCK, GRID>>>(dcolcand,dg1,dg2,dcols,dcnt,dsymbols);
        cudaDeviceSynchronize();        
	checkCudaErrors(cudaGetLastError());
}

void d_update_colcand(bool *dcolcand, discrete *dg1, discrete *dg2, int dcols)
{
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        update_colcand<<<BLOCK, GRID>>>(dcolcand, dg1, dg2, dcols);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
}

void d_intersect_row(bool *dcolcand, discrete *dg1, discrete *dg2, int dcols, int dcnt)
{        dim3 BLOCK(BLOCKsize,1,1);
         dim3 GRID(rows/BLOCKsize, 1,1);
         intersect_row<<<BLOCK, GRID>>>(dcolcand,dg1,dg2,dcols,dcnt);
        cudaDeviceSynchronize();
         checkCudaErrors(cudaGetLastError());
}
void d_get_pvalue(continuous da, int db, long double dpval){
        dim3 BLOCK(BLOCKsize,1,1);
        dim3 GRID(rows/BLOCKsize, 1,1);
        get_pvalue<<<BLOCK, GRID>>>(da, db, dpval);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
}

void d_clusprofile_init(int dcols, int dsigma, int *dprofile)
{
	dim3 BLOCK(BLOCKsize, 1,1);
	dim3 GRID(rows/BLOCKsize,1,1);
	clusprofile_init<<<BLOCK,GRID>>>(dcols,dsigma,dprofile);
	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());
}
