#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <htslib/sam.h>
#include <time.h>
#include "mt19937ar.h"
#include "filter.h"

size_t *GetSeqInt(size_t start, size_t end, size_t step);

size_t *SampleInt(size_t *arrayIn, size_t nTotal, size_t nSample, unsigned int replace, unsigned int seed);

void PrintArrayInt(size_t *array, size_t n);

int vsI(const void *a, const void *b);

node *cell_counts(gzFile R1_file, size_t len_cellbarcode, size_t len_umi);

typedef struct CB_node {
    node *CR;
    char *CB;
    struct CB_node *left;
    struct CB_node *right;
} CB_node;

CB_node *insert_CB_node(CB_node *root, char *CB, char *CR);
void free_CB_node(CB_node *root);
void print_CB_node(CB_node *root, gzFile fp);
CB_node *read_bam(char *bam_file);
void extract_bam(char *bam_file, const char *tag, int type);
#endif