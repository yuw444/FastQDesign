#define _POSIX_C_SOURCE 200809L

#include "unity.h"
#include "../../src/utils.h"
#include "../../src/filter.h"
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#define DATA_PATH "../../data/"

void setUp(void) {}
void tearDown(void) {}

void test_GetSeqInt_positive_step(void)
{
    size_t *seq = GetSeqInt(0, 10, 2);
    TEST_ASSERT_NOT_NULL(seq);
    TEST_ASSERT_EQUAL(0, seq[0]);
    TEST_ASSERT_EQUAL(2, seq[1]);
    TEST_ASSERT_EQUAL(4, seq[2]);
    TEST_ASSERT_EQUAL(6, seq[3]);
    TEST_ASSERT_EQUAL(8, seq[4]);
    TEST_ASSERT_EQUAL(10, seq[5]);
    free(seq);
}

void test_GetSeqInt_single_element(void)
{
    size_t *seq = GetSeqInt(5, 5, 1);
    TEST_ASSERT_NOT_NULL(seq);
    TEST_ASSERT_EQUAL(5, seq[0]);
    free(seq);
}

void test_SampleInt_without_replacement(void)
{
    size_t arrayIn[] = {1, 2, 3, 4, 5};
    size_t *sample = SampleInt(arrayIn, 5, 3, 0, 42);
    TEST_ASSERT_NOT_NULL(sample);
    for (int i = 0; i < 3; i++) {
        TEST_ASSERT_TRUE(sample[i] >= 1 && sample[i] <= 5);
    }
    free(sample);
}

void test_SampleInt_with_replacement(void)
{
    size_t arrayIn[] = {1, 2, 3};
    size_t *sample = SampleInt(arrayIn, 3, 10, 1, 42);
    TEST_ASSERT_NOT_NULL(sample);
    for (int i = 0; i < 10; i++) {
        TEST_ASSERT_TRUE(sample[i] >= 1 && sample[i] <= 3);
    }
    free(sample);
}

void test_vsI_comparison(void)
{
    size_t a = 5;
    size_t b = 10;
    TEST_ASSERT_TRUE(vsI(&a, &b) < 0);
    TEST_ASSERT_TRUE(vsI(&b, &a) > 0);
    TEST_ASSERT_EQUAL(0, vsI(&a, &a));
}

void test_insert_CB_node(void)
{
    CB_node *root = NULL;
    root = insert_CB_node(root, "CB1", "CR1");
    root = insert_CB_node(root, "CB2", "CR2");
    root = insert_CB_node(root, "CB1", "CR3");
    
    TEST_ASSERT_NOT_NULL(root);
    TEST_ASSERT_EQUAL_STRING("CB1", root->CB);
    TEST_ASSERT_NOT_NULL(root->CR);
    
    free_CB_node(root);
}

void test_free_CB_node(void)
{
    CB_node *root = NULL;
    root = insert_CB_node(root, "CB1", "CR1");
    root = insert_CB_node(root, "CB2", "CR2");
    
    free_CB_node(root);
}

void test_print_CB_node(void)
{
    CB_node *root = NULL;
    root = insert_CB_node(root, "CB1", "CR1");
    
    gzFile fp = gzopen(DATA_PATH "test_output.txt.gz", "wb");
    TEST_ASSERT_NOT_NULL(fp);
    
    print_CB_node(root, fp);
    gzclose(fp);
    
    free_CB_node(root);
}

void test_cell_counts_with_fastq(void)
{
    gzFile fp = gzopen(DATA_PATH "test_R1.fastq.gz", "rb");
    if (!fp) {
        TEST_FAIL();
        return;
    }
    
    node *root = cell_counts(fp, 16, 12);
    gzclose(fp);
    
    if (root) {
        TEST_ASSERT_TRUE(root->count >= 0);
        free_tree_node(root);
    }
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_GetSeqInt_positive_step);
    RUN_TEST(test_GetSeqInt_single_element);
    RUN_TEST(test_SampleInt_without_replacement);
    RUN_TEST(test_SampleInt_with_replacement);
    RUN_TEST(test_vsI_comparison);
    RUN_TEST(test_insert_CB_node);
    RUN_TEST(test_free_CB_node);
    RUN_TEST(test_print_CB_node);
    RUN_TEST(test_cell_counts_with_fastq);
    return UNITY_END();
}
