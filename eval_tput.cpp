#include <bits/stdc++.h>
#include <random>
#include <unistd.h>
#include "hashutil.h"
#include "cuckoo.h"

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a)  cerr << #a << " = " << a << endl

using namespace std;

const int numKeys = 20000000; // number of prepared insert/lookup keys
const int maxSteps = 200; // threshold for kick steps

vector<int> insKey, lupKey;
mt19937 rd;

template <typename fp_t, int fp_len>
void evaluate(Filter<fp_t, fp_len> &filter,
			  const char *filter_name)
{
	double insTime = 0, tput = 0, fac = 0;
	const int runTimes = 10;
	
	for (int run = 0; run < runTimes; ++run) {
		filter.clear();
		lupKey.clear(), insKey.clear();
		for (int i = 1; i <= numKeys; ++i)
			insKey.push_back(rd());
		for (int i = 1; i <= numKeys; ++i)
			lupKey.push_back(rd());

		int insStart = clock();
		for (int key : insKey) {
			if (filter.insert(key) == 1)
				break;
		}
		int insEnd = clock();
		insTime += insEnd - insStart;

		int lupStart = clock(), curTP = 0;
		for (int key : lupKey) {
			if (curTP % 10000 == 0 && clock() - lupStart >= CLOCKS_PER_SEC)
				break;
			filter.lookup(key);
			++curTP;
		}
		fac += filter.get_load_factor();
		tput += curTP;
	}

	int m = filter.n, b = filter.m;
	printf("Testing %s: m = %d, b = %d\n", filter_name, m, b);
	printf("    construction time = %.5f\n", insTime / CLOCKS_PER_SEC / runTimes);
	printf("    lookup throughput = %.0f\n", tput / runTimes);
	printf("    load factor = %.3f\n", fac / runTimes);
}

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed  = 20180804;
	int cmd_m = 18;
	while ((c = getopt(argc, argv, "r:l:")) != EOF) {
		switch (c) {
		case 'r':
			seed = atoi(optarg);
			break;
		case 'l':
			cmd_m = atoi(optarg);
			break;
		}	
	}
	rd.seed(seed);

	int m = 1 << cmd_m, b = 4; // #bucket, #entry per bucket
	CuckooFilter<uint8_t, 8> xor_filter;
	xor_filter.init(m, b, maxSteps);
	evaluate(xor_filter, "XOR-CuckooFilter");

	MyFilter<uint8_t, 8> add_filter;
	add_filter.init(m, b, maxSteps);
	evaluate(add_filter, "ADD-CuckooFilter");

	return 0;	
}
