#define _POSIX_C_SOURCE 200809L

#include "utils.h"
#include "filter.h"
#include <time.h>

size_t *GetSeqInt(size_t start, size_t end, size_t step)
{
  if (step == 0)
  {
    printf("Step cannot be 0");
    exit(1);
  }

  if ((end - start) * step < 0)
  {
    printf("Step must be positive when start < end, or negative when start > end");
    exit(1);
  }

  size_t *seq = (size_t *)calloc((size_t)(end + step -1 - start)/step + 1, sizeof(size_t));

  size_t i = 0;
  for (size_t j = start; j <= end; j += step)
  {
    seq[i] = j;
    i++;
  }

  return seq;
}

size_t *SampleInt(size_t *arrayIn, size_t nTotal, size_t nSample, unsigned int replace, unsigned int seed)
{

  init_genrand(seed);

  size_t *sampleOut = calloc(nSample, sizeof(size_t));

  if (replace == 0)
  {
    if (nSample > nTotal)
    {
      free(sampleOut);
      printf("Sample size must be smaller than population size when sampling without replacement.");
      exit(1);
    }

    size_t *arrayInCopy = calloc(nTotal, sizeof(size_t));
    memcpy(arrayInCopy, arrayIn, nTotal * sizeof(size_t));

    if (nTotal == nSample) {
        free(sampleOut);
        return arrayInCopy;
    }

    for (size_t i = 0; i < nSample; i++)
    {
      size_t index = genrand_int32() % nTotal;
      sampleOut[i] = arrayInCopy[index];
      if (index != (nTotal - 1))
      {
        arrayInCopy[index] = arrayInCopy[nTotal - 1];
      }
      nTotal--;
    }

    free(arrayInCopy);
  }
  else
  {
    for (size_t i = 0; i < nSample; i++)
    {
      sampleOut[i] = arrayIn[genrand_int32() % nTotal];
    }
  }

  return sampleOut;
}

void PrintArrayInt(size_t *array, size_t n)
{
  for (size_t i = 0; i < n; i++)
  {
    printf("%d ", array[i]);
  }
  printf("\n");
}

int vsI(const void *a, const void *b)
{
  return (*(size_t *)a - *(size_t *)b);
}

node *cell_counts(gzFile R1_file, size_t len_cellbarcode, size_t len_umi)
{
    node *root = NULL;

    while (1)
    {
        fastq *R1_block = get_fastq(R1_file);
        if (R1_block == NULL)
        {
            break;
        }
        char *cell_barcode_UMI = substring(R1_block->seq, 0, len_cellbarcode + len_umi);
        root = insert_tree(root, cell_barcode_UMI);
        free(cell_barcode_UMI);
        free_fastq(R1_block);
    }

    return root;
}

CB_node *insert_CB_node(CB_node *root, char *CB, char *CR)
{
    if (root == NULL)
    {
        root = malloc(sizeof(CB_node));
        root->CB = strdup(CB);
        root->CR = new_node(CR);
        root->left = NULL;
        root->right = NULL;
        return root;
    }
    int cmp = strcmp(CB, root->CB);
    if (cmp < 0)
    {
        root->left = insert_CB_node(root->left, CB, CR);
    }
    else if (cmp > 0)
    {
        root->right = insert_CB_node(root->right, CB, CR);
    }
    else
    {
        root->CR = insert_tree(root->CR, CR);
    }

    return root;
}

void free_CB_node(CB_node *root)
{
    if (root == NULL)
    {
        return;
    }

    free(root->CB);
    free_tree_node(root->CR);

    free_CB_node(root->left);
    free_CB_node(root->right);
    free(root);
}

void print_CB_node(CB_node *root, gzFile fp)
{
    if (root == NULL)
    {
        return;
    }

    gzprintf(fp, "%s;", root->CB);
    print_tree_same_row(root->CR, fp);
    gzprintf(fp, "\n");

    print_CB_node(root->left, fp);
    print_CB_node(root->right, fp);
}

CB_node *read_bam(char *bam_file)
{
    samFile *bam_reader = hts_open(bam_file, "r");

    if (bam_reader == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open bam file %s\n", bam_file);
        exit(1);
    }

    bam_hdr_t *bam_header = sam_hdr_read(bam_reader);
    bam1_t *bam_record = bam_init1();

    CB_node *root = NULL;
    long unsigned int read_count = 0;

    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        uint8_t *cb = bam_aux_get(bam_record, "CB");
        uint8_t *cr = bam_aux_get(bam_record, "CR");

        if (cb != NULL)
        {
            char *CB = (char *)malloc(18 + 1);
            char *CR = (char *)malloc(16 + 1);
            strcpy(CB, bam_aux2Z(cb));
            strcpy(CR, bam_aux2Z(cr));
            root = insert_CB_node(root, CB, CR);
            free(CB);
            free(CR);
        }

        read_count++;
        if (read_count % 10000000 == 0)
        {
            time_t t = time(NULL);
            struct tm *tm = localtime(&t);
            char s[64];
            strftime(s, sizeof(s), "%c", tm);
            printf("%s: ", s);
            printf("Processed %lu reads \n", read_count);
        }
    }

    printf("Processed all %lu reads\n", read_count);
    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_reader);

    return root;
}

void extract_bam(char *bam_file, const char *tag, int type)
{
    samFile *bam_reader = hts_open(bam_file, "r");

    if (bam_reader == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open bam file %s\n", bam_file);
        exit(1);
    }

    bam_hdr_t *bam_header = sam_hdr_read(bam_reader);
    bam1_t *bam_record = bam_init1();

    node *root = NULL;
    long total_count = 0;
    long valid_count = 0;

    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        total_count++;
        if (total_count % 10000000 == 0)
        {
            time_t t = time(NULL);
            struct tm *tm = localtime(&t);
            char s[64];
            strftime(s, sizeof(s), "%c", tm);
            printf("%s: ", s);
            printf("Processed %lu reads \n", total_count);
        }

        uint8_t *tag_ptr = bam_aux_get(bam_record, tag);
        char *tag_str = NULL;
        if (type){
            tag_str = calloc(1, sizeof(int));
        }

        if (tag_ptr != NULL)
        {
            valid_count++;
            switch (type)
            {
            case 0:
                root = insert_tree(root, bam_aux2Z(tag_ptr));
                break;
            case 1:
                sprintf(tag_str, "%d", bam_aux2i(tag_ptr));
                root = insert_tree(root, tag_str);
                free(tag_str);
                break;
            }
        }
    }

    FILE *fp = fopen("tag_summary.csv", "w");
    print_tree(root, fp);
    fclose(fp);

    free_tree_node(root);

    printf("Processed all %lu reads\n", total_count);
    printf("Valid reads: %lu\n", valid_count);
    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_reader);
}