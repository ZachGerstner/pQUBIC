#include "cluster.h"
#include <cuda.h>
#include <cuda-runtime.h>

static int compare_int(const void *a, const void *b)
{
	const int *da = a;
	const int *db = b;
	return (*da < *db)?-1:(*da != *db);
}

__device__ void update_colcand(bool *dcolcand, discrete *dg1, discrete *dg2, int dcols)
{
	int tx = threadIdx.x;x`
	for(int i = 0; i < tx; ++i)
		if((i < dcols) && dcolcand[i] && (dg1[i] != dg2[i]))
			dcolcand[i] = FALSE; 
}

__global__ void intersect_row(const bool *dcolcand, discrete *dg1, discrete *dg2, int dcols)
{
	int tx = threadIdx.x;
	int cnt = 0;
	for(int i = 0; i < tx; ++i)
	{
		if((i < dcols) && dcolcand[i] && (dg1[i] == dg2[i]) && (dg1[i] != 0))
			cnt++;
	}
	__syncthreads();
	return cnt;
}

__global__ static void reverse_row(const bool *dcolcand, discrete *dg1, discrete *dg2, int dcols, int dcnt)
{
	int tx = threadIdx.x;
	int cnt = 0;
	for(int i = 0; i < tx; ++i)
	{
		if((i < dcols) && dcolcand[i] && (symbols[dg1[i]] == -symbols[dg2[i]]))
			cnt++;
	}
	__syncthreads();
	dcnt = cnt;
}

__global__ static void seed_current_modify(const discrete *ds, bool *dcolcand, int *dcnt,int dthreshold,  int dcomponents, int dcols, int dsigma,int dtolerance, int *dprofile)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	int indy = blockDim.y * blockIdx.y + threadIdx.y;
	int k, flag, n;
	int threashold = ceil(dcomponents * dtolerance);
	discrete ss;
	*dcnt = 0;
	for(int i = 0; i < indx; ++i)
	{
		if(i < dcols)
		{
			flag = 0; ss = s[i];
			for(k = 1; k < indy; ++k)
			{
				if(k < dsigma)
				{
					n = dprofile[i][k];
					if ( k == ss)
						++n;
					if(n >= dthreshold)
					{
						flag = k;
						break;
					}
				}	
			}
			__syncthreads();
			if(flag)
			{
				*(cnt)++;
				dcolcand[i] = TRUE;
			} 
		}
	}
}

__global__ void bi_colcand_init(bool *dcolcand, int dcols, discrete *dg1, discrete *dg2)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	for(int i = 0; i < indx; ++i)
		if(i < dcols)
			dcolcand[i] = FALSE;
	for(int i = 0; i < indx; ++i)
		if(i < dcols && (dg1[i] == dg2[i]) && (dg1[i] != 0))
			dcolcand[i] = TRUE;
} 

__global__ void bi_candidate_100row(bool *dcandidates, int dtop, int *darr_rows)
{
	int indx = blockDim.x * blockIdx.x + threadIdx.x;
	for(int i = 0; i < indx; ++i)
	{
		if(darr_rows[i] < dtop)
			dcandidates[i] = FALSE;
	}
}

__global__ void d_profile_increment(Block *dbb, int db, int *dprofiles)
{
	int indx = blockIdx.x * blockDim.x + threadIdx.x;
        for ( i = 0; i < index; i++)
	{
		if(i < dbb[db]->block_rows);
        		dprofiles[dsItem(dbb[dbb]->genes,i)]++;
	}
}

static bool check_seed(Edge *e, Block **bb, const int block_id)
/*check whether current edge can be treat as a seed*/
{
        int profiles[rows];
        int dprofiles[rows];
	int i,b1,b2,b3,db1,db2;
	Block **dbb;
        bool fg = FALSE;
        b1 = b2 = -1;

	//CudaMalloc
	checkCudaErrors(cudaMalloc(dbb, sizeof(Block)*rows));
	checkCudaErrors(cudaMalloc(db1, sizeof(int)));
	checkCudaErrors(cudaMalloc(db2, sizeof(int)));
	checkCudaErrors(cudaMalloc(dprofiles, sizeof(int)*rows));

        for (i = 0; i < block_id; i++)
                if ( isInStack(bb[i]->genes,e->gene_one) && isInStack(bb[i]->genes, e->gene_two) )
                        return FALSE;

        for ( i = 0; i < rows; i++) profiles[i] = 0;
        fg = FALSE;
        for ( i = 0; i < block_id; i++)
                if ( isInStack(bb[i]->genes, e->gene_one) )
                {
                        fg = TRUE;
                        break;
                }
        if (fg)
                b1 = i;
        fg = FALSE;
        for ( i = 0; i < block_id; i++)
                if ( isInStack(bb[i]->genes, e->gene_two) )
                {
                        fg = TRUE;
                        break;
                }
        if (fg)
                b2 = i;
        if ( (b1 == -1)||(b2 == -1) )
                return TRUE;
        else
        {
		//moving values to gpu and calling gpu functions
		checkCudaErrors(cudaMemcpy(dbb, &bb, sizeof(BLOCK)*rows, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(db1, b1, sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(db2, b2, sizeof(int), cudaMemcpyHostToDevice));
		d_profile_wrap(dbb, db2, db1, dprofiles);
		//returning values from gpu
		checkCudaErrors(cudaMemcpy(&profiles, dprofiles, sizeof(int)*rows, cudaMemcpyDeviceToHost));
		for(i = 0; i < rows; ++i)
			if(profiles[i] > 1)
				return FALSE;
                b3 = MAX(bb[b1]->block_cols, bb[b2]->block_cols);
                if ( e->score <b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/ )
                        return FALSE;
                else
                        return TRUE;
        }
        err("never see this message\n");
        return FALSE;
}


void block_init(Edge *e, Block *b, struct dyStack *genes, struct dyStack *scores, bool *candidates,
				const int cand_threshold, int *components, struct dyStack *allincluster, 
				long double *pvalues)
{
	int i,score,dcols,top,dtop, cnt = 0,dnt, cnt_all = 0, pid = 0;
	continuous cnt_ave = 0, row_all = rows;
	long double pvalue;
	int max_cnt, max_i;
	int *arr_rows, *arr_rows_b;
	AllocArray(arr_rows, rows);
	AllocArray(arr_rows_b, rows);
	bool *colcand, *dcolcand;
	AllocArray(colcand, cols)
	discrete *g1, *g2, *dg2, *dg1;
	g1 = arr_c[dsItem(genes, 0)];
	g2 = arr_c[dsItem(genes, 1)];
	
	//Cuda Mallocs
	checkCudaErrors(cudaMalloc());

	//cuda Memcpy calls for initializing dcolcand
	cudaCheckErrors(cudaMemcpy(dcolcand, colcand, sizeof(bool)));
	cudaCheckErrors(cudaMemcpy(dg1, g1, sizeof(discrete)));
	cudaCheckErrors(cudaMemcpy(dg2, g2, sizeof(discrete)));
	cudaCheckErrors(cudaMemcpy(dcols, cols, sizeof(int)));
	//Wrapper function for colcand init
	bi_colcand_wrap(dcolcand, dg1, dg2, dcols);
	//cuda Memcpy back to host
	cudaCheckErrors(cudaMemcpy(colcand, dcolcand));
	for(i = 0; i < rows; ++i)
	{
		arr_rows[i] = intersect_row_wrap(dcolcand,dg1, darr_c[i], dcols);
		arr_rows_b[i] = arr_rows[i];
	}
	__syncthreads();
	 /*we just get the largest 100 rows when we initial a bicluster because we believe that 
         * the 100 rows can characterize the structure of the bicluster 
         * btw, it can reduce the time complexity*/
        if (rows > 100)
        {
                qsort(arr_rows_b, rows, sizeof *arr_rows, compare_int);
                top = arr_rows_b[rows -100];
		//wrapper memcpies and call
                cudaCheckErrors(cudaMemcpy(dtop, top, sizeof(int)));
		cand_false_wrap(drows, darr_rows, dcandidates, dtop);
		checkCudaErrors(cudaMemcpy(candidates, dcandidates, sizeof(bool)));
        }	
	//calculate the condition low bound for the current seed
	int cutoff = floor(0.05 * drows);
	db->cond_low_bound = arr_rows_b[rows-cutoff];

	while(*components < rows)
	{
		max_cnt = -1;
		max_i = -1;
		(*dcomponents)++;
		cnt_all = 0;
		cnt_ave = 0;

		//function of controlling the bicluster pval
		for(i = 0; i < rows; ++i)
		{
			if(!candidates[i]) continue;
			if(po->IS_list && !sublist[i]) continue;
			intersect_row_wrap(dcolcand, dg1, darr_c[i], dcols, dcnt); 
			checkCudaErrors(cudaMemcpy(cnt,dcnt, sizeof(int), cudaMemcpyDeviceToHost));
			cnt_all +=cnt;
			if(cnt< dcand_threshold)
				dcandidates = FALSE;
			if(cnt > max_cnt)
			{
				max_cnt = cnt;
				max_i = i;
			}
		}
		cnt_ave = cnt_all/row_all;
		pvalue = get_pvalue<<<BLOCKsize,GRIDsize>>>(cnt_ave, max_cnt);
		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());
		if(dpo->IS_cond)
		{
			if(max_cnt < dpo->COL_WIDTH || max_i < 0 || max_cnt < db->cond_low_bound)break;		
		}
		else
			if(max_cnt < po->COL_WIDTH || max_i < 0) break;
		if(dpo->IS_area)
			score = *dcomponents * max_cnt;
		else
			score = MIN(*dcomponent, max_cnt);
		if(score > db->score)
			db->sore = score;
		if(pvalue < db->pvalue)
			db->pvalue = pvalue;
		dsPush(genes, max_i);
		dsPush(scores, score);
		pvalues[pid++] = pvalue;
		update_colcand
	}
}

int cluster(FILE *fw, Edge **el, int n)
{
	int block_id = 0;
	Block **bb;
	int allocated = po->SCH_BLOCK;
	AllocArray(bb, allocated);
	Edge *e;
	Block *b;
	struct dyStack *genes, *scores, *b_genes, *allincluster;
	int i,j,k,components;
	AllocArray(profile, cols);
	for(j = 0l j < cols; ++j)
	{
		allocArray(profile[j], sigma);
	}
	genes = dsNew(rows);
	scores = dsNew(rows);
	allincluster = dsNew(rows);

	long double *pvalues;
	AllocArray(pvalues, rows);
	bool *candidates;
	AllocArray(candidates, rows);
	
	e = *el;
	i = 0;
	while (++i <n)
	{
		e = *el++;
		/*check if both genes already enumerated in previous blocks*/
		bool flag = TRUE;
		/*speed up the program if the rows abigger than 200*/
		if(rows > 250)
		{
			  if ( isInStack(allincluster,e->gene_one) && isInStack(allincluster,e->gene_two) )
                                flag = FALSE;
                        else if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
                                flag = FALSE;
                        else if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
                                flag =FALSE;
		}
		 else
                {
                        flag = check_seed(e, bb, block_id);
                        if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
                                flag = FALSE;
                        if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
                                flag = FALSE;
                }
                if (!flag) continue;
		d_profile_init_wrap(dprofile, dcols, dsigma);

		/*you must allocate a struct if you want to use the pointers related to it*/
		AllocVar(b);
		/*initial the b->score*/
		b->score = MIN(2, e->score);
		/*initial the b->pvalue*/
		b->pvalue = 1;

		/*initializing the stacks genes and scores*/
		int ii;
		dsClear(genes);
		dsClear(scores);
		for(ii = 0; ii < rows; ++ii)
		{
			dsPush(genes, -1);
			dsPush(scores, -1);
		}	
		dsClear(genes);
		dsClear(scores);
		
		dsPush(genes, e->gene_one);
		dsPush(genes, e->gene_two);
		dsPush(scores, 1);
		dsPush(scores, b->score);

		/*branch-and-cut condition for seed expandition*/
		int cand_threshold = floor(po->COL_WIDTH * po->TOLERANCE);
		if(cand_threshold < 2)
			cand_threshold = 2;

		/*maintain a candidate list to avoid looping through all rows*/
		//copy, call, and retrieve cuda values and functions
		checkCudaErrors(cudaMemcpy(dcandidates, candidates, sizeof(), cudaMemcpyHostToDevice));
		d_cand_init_wrap(drows, dcandidates);
		checkCudaErrors(cudaMemcpy(candidates, dcandidates, sizeof(), cudaMemcpyDeviceToHost));
		candidates[e->gene_one] = condidates[e->gene_two] = FALSE;
		components = 2;
		
		/* expansion step, generate a bicluster without noise */
                block_init(e, b, genes, scores, candidates, cand_threshold, &components, allincluster, pvalues);

                /* track back to find the genes by which we get the best score*/
                for(k = 0; k < components; k++)
                {
                        if (po->IS_pvalue)
                                if ((pvalues[k] == b->pvalue) &&(k >= 2) &&(dsItem(scores,k)!=dsItem(scores,k+1))) break;
                        if ((dsItem(scores,k) == b->score)&&(dsItem(scores,k+1)!= b->score)) break;
                }
                components = k + 1;
	
		//copy, retrieve, and call cuda values and functions
		checkCudaErrors(cudaMemcpy(drows, rows, sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(dcandidates, candidates, sizeof(), cudaMemcpyHostToDevice));
		d_cand_init_wrap(drows, dcandidates);
		checkCudaErrors(cudaMemcpy(candidates, dcandidates, sizeof(), cudaMemcpyDeviceToHost));

		/*add columns satisfy the conservative r*/
		//cuda stuff
		d_seed_wrap(darr_c[dsItem(genes, k)], dcolcand, &dcnt; dcomponents);
		
		/*add some new possible genes*/
		int m_cnt;
		
		
	}

}
