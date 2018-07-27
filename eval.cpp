#include <bits/stdc++.h>
#include <random>
#include <unistd.h>
#include "hashutil.h"
#include "filter.h"

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a)  cerr << #a << " = " << a << endl

using namespace std;

void evaluate(Filter &filter,
			  const char *filter_name,
			  const vector<int> &insKey,
			  const vector<int> &lupKey)
{
	int fail_insert = 0;
	int false_positive = 0;
	int false_negative = 0;
	
	for (int key : insKey)
		if (filter.insert(key) == 1) // insertion failed
			++fail_insert;

	for (int key : insKey)
		if (filter.lookup(key) == false) // false negative
			++false_negative;

	for (int key : lupKey)
		if (filter.lookup(key) == true) // false positive
			++false_positive;

	int n = insKey.size(), q = lupKey.size(), m = filter.n, b = filter.m;
	printf("testing %s: n = %d, q = %d, m = %d, b = %d\n",
		   filter_name, n, q, m, b);
    printf("fail insertion rate = %.2f%%\n", fail_insert * 100.0 / n);
	printf("false negative rate = %.2f%%\n", false_negative * 100.0 / n);
	printf("false positive rate = %.2f%%\n", false_positive * 100.0 / q);
	printf("load factor = %.2f\n", filter.get_load_factor());
}

const int maxSteps = 200; // threshold for kick steps
const int lookupNum = 1000000; // loopup operation numbers
vector<int> insKey, lupKey;

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed  = 0;
	int cmd_n = 100000;
	FILE *fp  = stdin;
	while ((c = getopt(argc, argv, "r:f:n:")) != EOF) {
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
				break;		
			default:
				break;
		}	
	}

	int n, m, b = 4; // #item, #bucket, #entry per bucket

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
		n = cmd_n;
		for (int i = 1; i <= n; ++i)
			insKey.push_back(rd());
		for (int i = 1; i <= lookupNum; ++i)
			lupKey.push_back(rd());
	}

	CuckooFilter xor_filter;
	m = 1 << (int)(ceil(log2(n / 4.0 / 0.95)));
	xor_filter.initialization(m, b, maxSteps);
	evaluate(xor_filter, "XOR-CuckooFilter", insKey, lupKey);

	MyFilter add_filter;
	m = (int) (n / 4.0 / 0.95);
	m += m & 1;
	add_filter.initialization(m, b, maxSteps);
	evaluate(add_filter, "ADD-CuckooFilter", insKey, lupKey);
	
}
