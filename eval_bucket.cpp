#include <bits/stdc++.h>
#include <random>
#include <unistd.h>
#include "hashutil.h"
#include "cuckoo.h"

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a)  cerr << #a << " = " << a << endl

using namespace std;

template <typename fp_t, int fp_len>
bool evaluate(Filter<fp_t, fp_len> &filter,
			  const char *filter_name,
			  const vector<int> &insKey,
			  const vector<int> &lupKey)
{
	int false_positive = 0;
	
	for (int key : insKey)
		if (filter.insert(key) == 1) // insertion failed
			return false;

	unordered_set<int> S(insKey.begin(), insKey.end());
	for (int key : lupKey)
		if (filter.lookup(key) == true && !S.count(key)) // false positive
			++false_positive;

	int n = insKey.size(), q = lupKey.size(), m = filter.n, b = filter.m;
	printf("Testing %s: n = %d, q = %d, m = %d, b = %d\n",
		   filter_name, n, q, m, b);
	printf("    false positive rate = %.3f%%\n", false_positive * 100.0 / q);
	printf("    load factor = %.4f\n", filter.get_load_factor());
	return true;
}

const int maxSteps = 200; // threshold for kick steps
const int numLookups = 1000000; // number of lookup operations
vector<int> insKey, lupKey;

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed = 0;
	int cmd_n = 100000;
	while ((c = getopt(argc, argv, "r:n:")) != EOF) {
		switch (c) {
		case 'r':
			seed = atoi(optarg);
			break;
		case 'n':
			cmd_n = atoi(optarg);
			break;
		default:
			break;
		}	
	}

	int n, q, m, b = 4; // #item, #lookup, #bucket, #entry per bucket
	
	mt19937 rd(seed);
	n = cmd_n;  q = numLookups;
	for (int i = 1; i <= n; ++i)
		insKey.push_back(rd());
	for (int i = 1; i <= q; ++i)
		lupKey.push_back(rd());

	m = 1 << (int)(ceil(log2(n / b)));
	for (int run = 0; ; ++run) {
		CuckooFilter<uint8_t, 8> xor_filter;
		xor_filter.init(m, b, maxSteps);
		if (evaluate(xor_filter, "XOR-CuckooFilter", insKey, lupKey))
			break;
		m <<= 1;
	}
	printf("------------------------------------------\n");

	int lb = n/b/2, rb = n/b;
	while (lb <= rb) {
		int mid = (lb + rb) / 2;
		MyFilter<uint8_t, 8> add_filter;
		add_filter.init(2*mid, b, maxSteps);
		if (evaluate(add_filter, "ADD-CuckooFilter", insKey, lupKey))
			rb = mid - 1;
		else
			lb = mid + 1;
	}

	return 0;	
}
