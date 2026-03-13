#define _POSIX_C_SOURCE 200809L

#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define DATA_PATH "/scratch/g/chlin/Yu/FastQDesign/data/"
#define BUILD_PATH "/scratch/g/chlin/Yu/FastQDesign/build/"

void setUp(void) {}
void tearDown(void) {}

void test_filter_command(void)
{
    system("mkdir -p " BUILD_PATH "test_filter_output");
    
    int ret = system("cd " BUILD_PATH "test_filter_output && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "../../fastF filter "
        "-I " DATA_PATH "test_I1.fastq.gz "
        "-R " DATA_PATH "test_R1.fastq.gz "
        "-r " DATA_PATH "test_R2.fastq.gz "
        "-w " DATA_PATH "whitelist.txt "
        "-l 16 -t 0.1 -s 42");
    
    TEST_ASSERT_TRUE(ret >= 0);
    
    system("rm -rf " BUILD_PATH "test_filter_output");
}

void test_freq_command(void)
{
    int ret = system("cd " BUILD_PATH " && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "./fastF freq "
        "-R " DATA_PATH "test_R1.fastq.gz "
        "-l 16 -w " DATA_PATH "whitelist.txt");
    
    TEST_ASSERT_TRUE(ret >= 0);
}

void test_extract_command(void)
{
    int ret = system("cd " BUILD_PATH " && "
        "LD_LIBRARY_PATH=/hpc/apps/htslib/1.22.1/lib "
        "./fastF extract "
        "-b " DATA_PATH "test.bam "
        "-t CB");
    
    TEST_ASSERT_EQUAL(0, ret);
    
    struct stat st;
    TEST_ASSERT_EQUAL(0, stat(BUILD_PATH "tag_summary.csv", &st));
    
    system("rm -f " BUILD_PATH "tag_summary.csv");
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_filter_command);
    RUN_TEST(test_freq_command);
    RUN_TEST(test_extract_command);
    return UNITY_END();
}
