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
	for (; ;) {
		int key = rd();
		if (filter.insert(key) == 1)
			break;
	}
	printf("Testing %s: m = %d, b = %d\n",
		   filter_name, filter.n, filter.m);
	printf("\tload factor = %.3f\n", filter.get_load_factor());
}

const int maxSteps = 200; // threshold for kick steps

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed  = 0;
	int cmd_m = 100000;
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

	int m = cmd_m, b = 4; // #item, #lookup, #bucket, #entry per bucket

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
