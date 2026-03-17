#define _POSIX_C_SOURCE 200809L

#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DATA_PATH "/scratch/g/chlin/Yu/FastQDesign/data/"
#define BUILD_PATH "/scratch/g/chlin/Yu/FastQDesign/build/"

void setUp(void) {}
void tearDown(void) {}

void test_filter_command(void)
{
    system("mkdir -p " BUILD_PATH "test_filter_output");
    
    int ret = system("cd " BUILD_PATH "test_filter_output && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "../fastF filter "
        "-I " DATA_PATH "test_I1.fastq.gz "
        "-R " DATA_PATH "test_R1.fastq.gz "
        "-r " DATA_PATH "test_R2.fastq.gz "
        "-w " DATA_PATH "whitelist.txt "
        "-l 16 -t 0.1 -s 42");
    
    TEST_ASSERT_EQUAL(0, ret);
}

void test_freq_command(void)
{
    system("mkdir -p " BUILD_PATH "freq_output");
    
    int ret = system("cd " BUILD_PATH " && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "./fastF freq "
        "-R " DATA_PATH "test_R1.fastq.gz "
        "-l 16 -o " BUILD_PATH "freq_output/");
    
    TEST_ASSERT_EQUAL(0, ret);
}

void test_extract_command(void)
{
    int ret = system("cd " BUILD_PATH " && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "./fastF extract "
        "-b " DATA_PATH "test.bam "
        "-t CB");
    
    TEST_ASSERT_EQUAL(0, ret);
}

void test_crb_command(void)
{
    int ret = system("cd " BUILD_PATH " && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "./fastF crb "
        "-b " DATA_PATH "test.bam "
        "-o " BUILD_PATH "crb_output.tsv.gz");
    
    TEST_ASSERT_EQUAL(0, ret);
}

void test_bam2db_command(void)
{
    system("mkdir -p " BUILD_PATH "bam2db_output");
    
    int ret = system("cd " BUILD_PATH " && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "./fastF bam2db "
        "-b " DATA_PATH "test.bam "
        "-f " DATA_PATH "features.tsv.gz "
        "-a " DATA_PATH "barcodes.tsv.gz "
        "-o " BUILD_PATH "bam2db_output/");
    
    TEST_ASSERT_EQUAL(0, ret);
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_filter_command);
    RUN_TEST(test_freq_command);
    RUN_TEST(test_extract_command);
    RUN_TEST(test_crb_command);
    RUN_TEST(test_bam2db_command);
    return UNITY_END();
}
