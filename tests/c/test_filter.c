#define _POSIX_C_SOURCE 200809L

#include "unity.h"
#include "../../src/filter.h"
#include <string.h>
#include <stdlib.h>

#define DATA_PATH "../../data/"

void setUp(void) {}
void tearDown(void) {}

void test_new_node(void)
{
    node *n = new_node("test");
    TEST_ASSERT_NOT_NULL(n);
    TEST_ASSERT_EQUAL_STRING("test", n->data);
    TEST_ASSERT_EQUAL(1, n->count);
    TEST_ASSERT_NULL(n->left);
    TEST_ASSERT_NULL(n->right);
    free(n->data);
    free(n);
}

void test_insert_tree(void)
{
    node *root = NULL;
    root = insert_tree(root, "apple");
    root = insert_tree(root, "banana");
    root = insert_tree(root, "cherry");
    
    TEST_ASSERT_NOT_NULL(root);
    TEST_ASSERT_EQUAL_STRING("apple", root->data);
    TEST_ASSERT_TRUE(root->count >= 1);
    
    free_tree_node(root);
}

void test_in_tree_found(void)
{
    node *root = NULL;
    root = insert_tree(root, "apple");
    root = insert_tree(root, "banana");
    root = insert_tree(root, "cherry");
    
    TEST_ASSERT_TRUE(in(root, "apple", 6));
    TEST_ASSERT_TRUE(in(root, "banana", 6));
    TEST_ASSERT_TRUE(in(root, "cherry", 6));
    
    free_tree_node(root);
}

void test_in_tree_not_found(void)
{
    node *root = NULL;
    root = insert_tree(root, "apple");
    root = insert_tree(root, "banana");
    
    TEST_ASSERT_FALSE(in(root, "grape", 6));
    
    free_tree_node(root);
}

void test_substring(void)
{
    char *str = "HelloWorld";
    char *sub = substring(str, 0, 5);
    TEST_ASSERT_EQUAL_STRING("Hello", sub);
    free(sub);
    
    sub = substring(str, 5, 5);
    TEST_ASSERT_EQUAL_STRING("World", sub);
    free(sub);
}

void test_construct_tree(void)
{
    char *whitelist[] = {"AAAA", "BBBB", "CCCC"};
    node *root = construct_tree(whitelist, 3);
    
    TEST_ASSERT_NOT_NULL(root);
    TEST_ASSERT_EQUAL_STRING("AAAA", root->data);
    
    free_tree_node(root);
}

void test_get_row(void)
{
    int rows = get_row(DATA_PATH "whitelist.txt");
    TEST_ASSERT_GREATER_THAN(0, rows);
}

void test_read_txt(void)
{
    char **lines = read_txt(DATA_PATH "whitelist.txt", 10);
    TEST_ASSERT_NOT_NULL(lines);
    TEST_ASSERT_TRUE(strlen(lines[0]) > 0);
    
    for (int i = 0; i < 10; i++) {
        free(lines[i]);
    }
    free(lines);
}

int main(void)
{
    UNITY_BEGIN();
    
    RUN_TEST(test_new_node);
    RUN_TEST(test_insert_tree);
    RUN_TEST(test_in_tree_found);
    RUN_TEST(test_in_tree_not_found);
    RUN_TEST(test_substring);
    RUN_TEST(test_construct_tree);
    RUN_TEST(test_get_row);
    RUN_TEST(test_read_txt);
    
    return UNITY_END();
}
