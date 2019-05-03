/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Jan. 25, 2010 
 * Usage: This is part of the bicluster package. Use, redistribute, modify
 *        without limitations.
 *
 * Produces two graphs sequentially, derived from microarray data.
 * 
 * The first graph is generated from the raw data where an edge is defined
 * as two genes having common components from same condition, the score on
 * the edge being the number of the same components. The edges in the first
 * graph, are used as vertices in the second graph, where edges are defined
 * as the common columns between two co-vertex edges in the first graph,
 * with scores defined as the number of same columns for all three genes.
 *
 */

#include "make_graph.h"
#include "utils.h"
#include <cuda.h>
#include <cuda-runtime.h>
/*we can reduce the HEAP_SIZE when the data contain so many genes so that memory is not enough*/

/**************************************************************************/

/* String intersection function without string copying, only numbers */
/*icaculate the weight of the edge in the first graph*/

__global__ static void str_instersect_r(const discrete *s1, const discrete *s2, int dcommon_cnt, int dcols)
{
	int tx = threadIdx.x;
	int common_cnt = 0, i;
	for(i = 0; i <tx; ++i)
	{
		if(i < dcols)
		{
			if(*s1 == *s2 && (*s1!=0))
				common_cnt++;
			s1++;
			s2++;
		}
	}
	__syncthreads();
	dcommon_cnt = common_cnt;
}

__global__ void seed_deduct(const discrete *s, int *dprofile, int dcols)
{
	int tx = threadIdx.x;
	int i;
	discrete ss;
	for(i = 0; i < tx; ++i)
	{
		if(i < dcols)
		{
			ss = s[i];
			dprofile[i][ss]--;
		}
	}

}

__global__ void seed_update(const discrete *s, int *dprofiles, int dcols)
{
	int tx = threadIdx.x;
	int i;
	for(i = 0; i < tx; i++)
	{
		if(i < dcols)
			dprofile[i][s[i]]++;	
	}
}


