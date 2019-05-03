#ifndef _CLUSTER_H
#define _CLUSTER_H

#include "struct.h"
#define VER 1.0
/* this structure holds the matching score for each row */
struct rowMatch
{
    int row_id;
    int matches;
};

/* struct */
extern void verboseDot();

/* make_graph */
extern void seed_update (const discrete *s);

/* write_block */
extern void scan_block (struct dyStack *gene_set, Block *b_ptr);
extern void print_bc (FILE *fw, Block *b, int num);

/* cluster.cu */
void bi_colcand_wrap(bool *dcolcand, int dcols, discrete *dg1, discrete *dg2);
void bi_candidate_wrap(bool *dcandidates, int dtop, int *darr_rows);
void d_profile_increment(Block **dbb, int *db1, int *db2, int *dprofiles);
void d_seed_wrap(const discrete *ds, bool *dcolcand, int *dcnt, int dcomponents, int dcols, int dsigma, double dtolerance, int **dprofile);
void d_reverse_wrap(const bool *dcolcand, discrete *dg1, discrete *dg2, int dcols, int dcnt);
void d_update_colcand(bool *dcolcand, discrete *dg1, discrete *dg2, int dcols);
void d_intersect_row(bool *dcolcand, discrete *dg1, discrete *dg2, int dcols, int dcnt);
void d_get_pvalue(continuous da, int db, long double dpval);
void d_profile_init(int *dprofiles, int *drows);
void d_colcand_init(bool *dcolcand, int dcols);

/* prototypes */
static int compare_int (const int *a, const int *b);
static void update_colcand(bool *colcand, discrete *g1, discrete *g2);
static int intersect_row(const bool *colcand, discrete *g1, discrete *g2);
static int reverse_row(const bool *colcand, discrete *g1, discrete *g2);
static void seed_current_modify (const discrete *s, bool *colcand, int* cnt, int components);
static bool check_seed(Edge *e, Block **bb, const int block_id);
static void print_params(FILE *fw);
int cluster (FILE *fw, Edge **el, int n);
static int report_blocks(FILE *fw, Block **bb, int num);
static void sort_block_list(Block **el, int n);
static void block_init(Edge *e, Block *b, struct dyStack *genes, struct dyStack *scores,
bool *candidates, const int cand_threshold,int *components, struct dyStack *allincluster, long double *pvalues);
long double get_pvalue (continuous a, int b);

bits16 **profile;

#endif
