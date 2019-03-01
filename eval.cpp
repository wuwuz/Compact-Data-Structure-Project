#include <bits/stdc++.h>
#include <random>
#include <unistd.h>
#include "hashutil.h"
#include "cuckoo.h"

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a)  cerr << #a << " = " << a << endl

using namespace std;

template <typename fp_t, int fp_len>
double evaluate(Filter<fp_t, fp_len> &filter,
			  const char *filter_name,
			  const vector<int> &insKey,
			  const vector<int> &lupKey)
{
    // return load factor
    //
	int fail_insert = 0;
	int false_positive = 0;
	int false_negative = 0;
	
	for (int key : insKey)
    {
        if (key == insKey[515])
            printf("now");
		if (filter.insert(key) == 1) // insertion failed
        {
            printf("insert false : %x\n", key);
			++fail_insert;
            //break;
        }
    }

	for (int key : insKey)
		if (filter.lookup(key) == false) // false negative
        {
            printf("lookup false : %x\n", key);
			++false_negative;
        }

	unordered_set<int> S(insKey.begin(), insKey.end());
	for (int key : lupKey)
		if (filter.lookup(key) == true && !S.count(key)) // false positive
			++false_positive;

	int n = insKey.size(), q = lupKey.size(), m = filter.n, b = filter.m;
	printf("testing %s: n = %d, q = %d, m = %d, b = %d memory = %d bytes\n",
		   filter_name, n, q, m, b, filter.memory_consumption);
    printf("fail insertion rate = %.3f%%\n", fail_insert * 100.0 / n);
	printf("false negative rate = %.3f%%\n", false_negative * 100.0 / n);
	printf("false positive rate = %.3f%%\n", false_positive * 100.0 / q);
	printf("load factor = %.3f\n", filter.get_load_factor());

    return filter.get_load_factor();
}

const int maxSteps = 500; // threshold for kick steps
vector<int> insKey, lupKey;

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed  = 0;
	int cmd_n = 100000;
	int cmd_q = 1000000;
	FILE *fp  = stdin;
	while ((c = getopt(argc, argv, "r:f:n:q:")) != EOF) {
		switch (c) {
		case 'r':
			seed = atoi(optarg);
			break;
		case 'f':
			fp = fopen(optarg, "r");
			seed = -1;
			break;
		case 'n':
			cmd_n = atoi(optarg);
            printf("cmd_n == %d\n", cmd_n);
			break;
		case 'q':
			cmd_q = atoi(optarg);
            printf("cmd_q == %d\n", cmd_q);
			break;
		default:
			break;
		}	
	}

	int n, q, m, b = 4; // #item, #lookup, #bucket, #entry per bucket

	// read testdata from file
	int key;
	if (seed == -1) {
		fscanf(fp, "%d", &n);
		for (int i = 1; i <= n; ++i) {
			fscanf(fp, "%d", &key);
			insKey.push_back(key);	
		}
		while (fscanf(fp, "%d", &key) != EOF)
			lupKey.push_back(key);
	}
	else {
		mt19937 rd(seed);
		n = cmd_n;  q = cmd_q;
		for (int i = 1; i <= n; ++i)
			insKey.push_back(rd());
		for (int i = 1; i <= q; ++i)
			lupKey.push_back(rd());
	}
    
    printf("testing with parameters : \n");
    printf("insert keys : %d\n", n);
    printf("lookup times: %d\n\n\n", q);

    /*
	StandardCuckooFilter<uint8_t, 8> standard_cuckoo_filter;
	m = 1 << (int)(ceil(log2(n / b / 0.95)));
	standard_cuckoo_filter.init(m, b, maxSteps);
	evaluate(standard_cuckoo_filter, "standard-CuckooFilter", insKey, lupKey);

	MyFilter<uint8_t, 8> add_filter;
	//m = 1 << (int)(ceil(log2(n / b / 0.95)));
	m = (int) (n / b / 0.95 + 1);
	m += m & 1;
	add_filter.init(m, b, maxSteps);
	evaluate(add_filter, "OLD-CuckooFilter", insKey, lupKey);
    
    */

    /*

    printf("maxSteps %d\n", maxSteps);
    double sum1 = 0;
    double sum2 = 0;
    int times = 50;

	mt19937 rd(seed);

    for (int i = 1; i <= times; i++)
    {
        printf("round : %d\n", i);
        insKey.clear();
        lupKey.clear();

		n = cmd_n;  q = cmd_q;
		for (int i = 1; i <= n; ++i)
			insKey.push_back(rd());
		for (int i = 1; i <= q; ++i)
			lupKey.push_back(rd());

        XorFilter<uint8_t, 8> xor_filter;
        //m = 1 << (int)(ceil(log2(n / b / 0.95)));
        m = n / b;
        //m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        xor_filter.init(m, b, maxSteps);
        double sb = evaluate(xor_filter, "xor-CuckooFilter", insKey, lupKey);
        printf("%.5f\n", sb);
        sum1 += sb;

        AddFilter<uint8_t, 8> new_filter;
        //m = 1 << (int)(ceil(log2(n / b / 0.95)));
        //m = (int) (n / b / 0.96 + 1);
        m = n / b;
        m += m & 1;
        new_filter.init(m, b, maxSteps);
        sum2 += evaluate(new_filter, "Morton-CuckooFilter", insKey, lupKey);
    }

    printf("%.5f %.5f\n", sum1 / times, sum2 / times);

    */

    //b = 16;

    /*
    XorFilter<uint16_t, 16> xor_filter;
    xor_filter.init(1000, b, maxSteps);
    xor_filter.test_bucket();

    XorFilter<uint8_t, 8> xor_filter_1;
    xor_filter_1.init(1000, b, maxSteps);
    xor_filter_1.test_bucket();
    */

    // uint16 : buggy...

    MortonAddFilter<uint8_t, 8> new_filter;
    //m = 1 << (int)(ceil(log2(n / b / 0.95)));
    //m = (int) (n / b / 0.95 + 1);
    //m = n / b;
    //m += m & 1;
    new_filter.init(n, b, maxSteps);
    evaluate(new_filter, "Morton-CuckooFilter", insKey, lupKey);

    XorFilter<uint16_t, 15> xor_filter;
    //m = 1 << (int)(ceil(log2(n / b / 0.95)));
    //m = n / b;
    m = (int) (n / b / 0.90 + 1);
    m += m & 1;
    xor_filter.init(m, b, maxSteps);
    evaluate(xor_filter, "xor-CuckooFilter", insKey, lupKey);
    //printf("%.5f\n", sb);
    //sum1 += sb;

    BloomFilter<uint8_t, 8> bloom_filter;
	m = (int) (n / b / 0.95 + 1);
	m += m & 1;
	bloom_filter.init(m, b, maxSteps);
	evaluate(bloom_filter, "bloom-CuckooFilter", insKey, lupKey);

	return 0;	
}
