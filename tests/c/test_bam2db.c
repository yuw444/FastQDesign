#define _POSIX_C_SOURCE 200809L

#include "unity.h"
#include "../../src/bam2db.h"
#include <string.h>
#include <stdlib.h>

void setUp(void) {}
void tearDown(void) {}

void test_encode_base(void)
{
    TEST_ASSERT_EQUAL(0b00, encode_base('A'));
    TEST_ASSERT_EQUAL(0b01, encode_base('C'));
    TEST_ASSERT_EQUAL(0b10, encode_base('G'));
    TEST_ASSERT_EQUAL(0b11, encode_base('T'));
    TEST_ASSERT_EQUAL(0xFF, encode_base('X'));
}

void test_encode_DNA(void)
{
    uint8_t *encoded = encode_DNA("ACGT");
    TEST_ASSERT_NOT_NULL(encoded);
    TEST_ASSERT_TRUE(encoded[0] != 0);
    free(encoded);

    encoded = encode_DNA("AAAA");
    TEST_ASSERT_NOT_NULL(encoded);
    free(encoded);

    encoded = encode_DNA("invalid");
    TEST_ASSERT_NULL(encoded);
}

void test_decode_DNA(void)
{
    uint8_t encoded[] = {0xFF, 0xFF};
    char *decoded = decode_DNA(encoded, 4);
    TEST_ASSERT_NOT_NULL(decoded);
    TEST_ASSERT_TRUE(strlen(decoded) == 4);
    free(decoded);
}

void test_encode_decode_roundtrip(void)
{
    char *original = "ACGTACGTACGT";
    uint8_t *encoded = encode_DNA(original);
    TEST_ASSERT_NOT_NULL(encoded);
    
    char *decoded = decode_DNA(encoded, strlen(original));
    TEST_ASSERT_NOT_NULL(decoded);
    TEST_ASSERT_EQUAL_STRING(original, decoded);
    
    free(encoded);
    free(decoded);
}

void test_hash_function(void)
{
    uint64_t h1 = hash("test", 4);
    uint64_t h2 = hash("test", 4);
    uint64_t h3 = hash("other", 5);
    
    TEST_ASSERT_EQUAL(h1, h2);
    TEST_ASSERT_TRUE(h1 != h3);
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_encode_base);
    RUN_TEST(test_encode_DNA);
    RUN_TEST(test_decode_DNA);
    RUN_TEST(test_encode_decode_roundtrip);
    RUN_TEST(test_hash_function);
    return UNITY_END();
}
