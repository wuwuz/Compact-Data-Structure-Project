#include <bits/stdc++.h>
#include <random>
#include <unistd.h>
#include "hashutil.h"
#include "cuckoo.h"

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a)  cerr << #a << " = " << a << endl

using namespace std;

mt19937 rd;

template <typename fp_t, int fp_len>
void evaluate(Filter<fp_t, fp_len> &filter,
			  const char *filter_name,
			  int numLookups)
{
	const int runTimes = 10;
	int false_positive = 0;
	for (int run = 0; run < runTimes; ++run) {
		unordered_set<int> S;
		for (; ;) {
			int key = rd(); S.insert(key);
			if (filter.insert(key) == 1)
				break;
		}

		for (int i = 0; i < numLookups; ++i) {
			int key = rd();
			if (filter.lookup(key) == true && !S.count(key))
				++false_positive;
		}
	}

	int q = numLookups, m = filter.n, b = filter.m;
	printf("Testing %s: m = %d, b = %d, q = %d\n",
		   filter_name, m, b, q);
	printf("    false positive rate = %.5f%%\n", false_positive * 100.0 / q / runTimes);
}

const int maxSteps = 200; // threshold for kick steps
vector<int> insKey, lupKey;

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed  = 37;
	int cmd_q = 1000000;
	while ((c = getopt(argc, argv, "r:q:")) != EOF) {
		switch (c) {
		case 'r':
			seed = atoi(optarg);
			break;
		case 'q':
			cmd_q = atoi(optarg);
			break;
		default:
			break;
		}	
	}
	rd.seed(seed);

	int q = cmd_q, m = 262144, b = 4; // #lookup, #bucket, #entry per bucket

	CuckooFilter<uint16_t, 16> xor_filter;
	xor_filter.init(m, b, maxSteps);
	evaluate(xor_filter, "XOR-CuckooFilter", q);

	MyFilter<uint16_t, 16> add_filter;
	add_filter.init(m, b, maxSteps);
	evaluate(add_filter, "ADD-CuckooFilter", q);

	return 0;	
}
