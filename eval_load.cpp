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
			  const char *filter_name)
{
	double aver_load_factor = 0;
	int runTimes = 25;
	for (int run = 0; run < runTimes; ++run) {
		for (; ;) {
			int key = rd();
			if (filter.insert(key) == 1)
				break;
		}
		aver_load_factor += filter.get_load_factor();
	}
	printf("Testing %s: m = %d, b = %d\n",
		   filter_name, filter.n, filter.m);
	printf("    load factor = %.3f\n", aver_load_factor / runTimes);
}

const int maxSteps = 200; // threshold for kick steps

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed  = 97;
	int cmd_m = 262144;
	while ((c = getopt(argc, argv, "r:m:")) != EOF) {
		switch (c) {
		case 'r':
			seed = atoi(optarg);
			break;
		case 'm':
			cmd_m = atoi(optarg);
			break;
		default:
			break;
		}	
	}
	rd.seed(seed);

	int m = cmd_m, b = 4; // #bucket, #entry per bucket

	if ((m & -m) == m) {
		CuckooFilter<uint8_t, 8> xor_filter;
		xor_filter.init(m, b, maxSteps);
		evaluate(xor_filter, "XOR-CuckooFilter");
	}

	MyFilter<uint8_t, 8> add_filter;
	add_filter.init(m, b, maxSteps);
	evaluate(add_filter, "ADD-CuckooFilter");

	return 0;	
}
