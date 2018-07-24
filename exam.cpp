#include <bits/stdc++.h>
#include "hashutil.h"
#include "filter.h"

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a) cerr << #a << " = " << a << endl

using namespace std;
const int N = 2000010;
const int K = 200; // threshold for block split

int exam[N];

int main()
{

    srand(12345);

    int max_item = 10000;
    int n = 1;
    for (; n * 4 <= max_item; ) n <<= 1;

    deln(n);
    int m = 4;

    CuckooFilter a;
    // MyFilter a;
    a.initialization(n, m, 200);

    for (int i = 1; i <= max_item * 2; i++) exam[i] = i;
    random_shuffle(exam + 1, exam + 1 + max_item * 2);

    int fail_insert = 0;
    int false_positive = 0;
    int false_negative = 0;

    for (int i = 1; i <= max_item; i++) 
    {
        if (a.insert(exam[i]) == 1) // insert fail
        {
            fail_insert++;
        }
    }

    // positive lookup

    for (int i = 1; i <= max_item; i++)
    {
        if (a.lookup(exam[i]) == false) // false neage
        {
            false_negative++;
        }
    }

    for (int i = max_item + 1; i <= max_item * 2; i++)
    {
        if (a.lookup(exam[i]) == true) // false positive
            false_positive++;
    }

    printf("testing : n = %d, m = %d\n", n, m);
    printf("insert item = %d\n", max_item);
    printf("insert fail rate = %.2f\n", fail_insert * 100.0 / max_item);
    printf("false negative = %.2f\n", false_negative * 100.0 / max_item);
    printf("false positive = %.2f\n", false_positive * 100.0 / max_item);
    printf("load factor = %.2f\n", a.get_load_factor());

    return 0;
}


