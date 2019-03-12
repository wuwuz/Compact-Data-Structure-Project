#include <bits/stdc++.h>
#include <random>
#include <unistd.h>
#include "hashutil.h"
#include "cuckoo.h"

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a)  cerr << #a << " = " << a << endl

using namespace std;

template <typename fp_t, int fp_len>
void evaluate(Filter<fp_t, fp_len> &filter,
			  const char *filter_name,
			  const vector<int> &insKey,
			  const vector<int> &lupKey,
              vector<double> &result,
              bool verbose = false)
{
    // return load factor
    //
	int fail_insert = 0;
	int false_positive = 0;
	int false_negative = 0;

    printf("insert\n");
	
	for (int key : insKey)
    {
		if (filter.insert(key) == 1) // insertion failed
        {
            //printf("insert false : %x\n", key);
			++fail_insert;
            break;
        }
        //printf("%.5f\n", filter.get_load_factor());
    }

    printf("positive lookup\n");
	for (int key : insKey)
    {
		if (filter.lookup(key) == false) // false negative
        {
            filter.lookup(key);
            //printf("lookup false : %x\n", key);
			++false_negative;
        }
    }

    printf("negative lookup\n");

	unordered_set<int> S(insKey.begin(), insKey.end());

	for (int key : lupKey)
		if (filter.lookup(key) == true && !S.count(key)) // false positive
			++false_positive;

    printf("negative lookup finish\n");

	int n = insKey.size(), q = lupKey.size(), m = filter.n, b = filter.m;

    if (verbose)
    {
        printf("testing %s: n = %d, q = %d, m = %d, b = %d memory = %d bytes\n",
               filter_name, n, q, m, b, filter.memory_consumption);
        printf("fail insertion rate = %.3f%%\n", fail_insert * 100.0 / n);
        printf("false negative rate = %.3f%%\n", false_negative * 100.0 / n);
        printf("false positive rate = %.3f%%\n", false_positive * 100.0 / q);
        printf("load factor = %.3f\n", filter.get_load_factor());
        printf("full bucket factor = %.3f\n", filter.get_full_bucket_factor());
    }

    //result : 
    //0 : memory usage
    //1 : fail insert rate
    //2 : false negative rate
    //3 : false positive rate
    //4 : load factor
    //5 : bits per item

    result.resize(10);
    result[0] = filter.memory_consumption;
    result[1] = fail_insert * 1.0 / n;
    result[2] = false_negative * 1.0 / n;
    result[3] = false_positive * 1.0 / q;
    result[4] = filter.get_load_factor();
    result[5] = filter.memory_consumption * 8.0 / n;
}

void test_memory_usage()
{
    int m;
    int b = 4;
    int maxSteps = 400;

    FILE *out = NULL;
    out = fopen("key-memory.csv", "w");

    fprintf(out, "key number, Semi-Sort-Cuckoo-Filter-12bit, Morton-Filter, Vacuum-Filter-12bit, Vacuum-Filter-11bit, Bloom-Filter-12bit\n");
    int max = 10000000;

    for (int n = 1000; n <= max; n += 1000)
    {
        if (n % (max / 100) == 0)
            printf("process : %d%\n", n / (max / 100));

        StandardCuckooFilter<uint16_t, 12> sscf;
        m = 1 << (int)(ceil(log2(n / b / 0.95)));
        sscf.init(m, b, maxSteps);

        MortonAddFilter<uint8_t, 8> mf;
        m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        mf.init(m, b, maxSteps);

        VacuumFilter<uint16_t, 12> vf;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        vf.init(m, b, maxSteps);

        VacuumFilter<uint16_t, 11> vf2;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        vf2.init(m, b, maxSteps);

        BloomFilter<uint16_t, 12> bf;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        bf.init(m, b, maxSteps);

        fprintf(out, "%d, %d, %d, %d, %d, %d\n", n, sscf.memory_consumption, mf.memory_consumption, vf.memory_consumption, vf2.memory_consumption, bf.memory_consumption);
    }

    fclose(out);
}

void test_all()
{
    FILE *out = NULL;
    out = fopen("all-performance.csv", "w");

    //result : 
    //0 : memory usage
    //1 : fail insert rate
    //2 : false negative rate
    //3 : false positive rate
    //4 : load factor
    //5 : bits per item

    fprintf(out, "key size, ");

    const int test_time = 11;
    const int key[test_time] = {(1 << 16) + 1, 10000, 100000, 200000, 300000, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000};
    const string name[4] = {"Vacuum", "Standard", "Morton", "Bloom"};
    const string metric[6] = {"memory", "fail insert", "false negative", "false positive", "load factor", "bits per item"};

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 4; j++)
        {
            string tmp = name[j] + " " + metric[i];
            if (i != 5 || j != 3)
                fprintf(out, "%s, ", tmp.c_str());
            else 
                fprintf(out, "%s\n", tmp.c_str());
        }

    for (int i = 0; i < test_time; i++)
    {
        int n = key[i];
        int q = 1000000;
        int maxSteps = 400;
        int b = 4;
        int seed = 1;
        int m;


        vector<int> insKey;
        vector<int> lupKey;

        mt19937 rd(seed);
        for (int i = 1; i <= n; ++i)
            insKey.push_back(rd());
        for (int i = 1; i <= q; ++i)
            lupKey.push_back(rd());

        vector<double> res[4];

        VacuumFilter<uint16_t, 11> vf;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        vf.init(m, b, maxSteps);
        evaluate(vf, name[0].c_str(), insKey, lupKey, res[0]);

        deln(m);

        StandardCuckooFilter<uint16_t, 11> cf;
        m = 1 << (int)(ceil(log2(n / b / 0.95)));
        cf.init(m, b, maxSteps);
        evaluate(cf, name[1].c_str(),  insKey, lupKey, res[1]);


        MortonAddFilter<uint8_t, 8> mf;
        m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        mf.init(m, b, maxSteps);
        evaluate(mf, name[2].c_str(),  insKey, lupKey, res[2]);

        printf("start\n");
        BloomFilter<uint16_t, 11> bf;
        m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        bf.init(m, b, maxSteps);
        printf("mid\n");
        evaluate(bf, name[3].c_str(),  insKey, lupKey, res[3]);
        printf("finish\n");

        fprintf(out, "%d, ", n);

        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 4; k++)
            {
                if (j != 5 || k != 3)
                    fprintf(out, "%.5f, ", res[k][j]);
                else
                    fprintf(out, "%.5f\n", res[k][j]);
            }
    }
}

const int maxSteps = 250; // threshold for kick steps
vector<int> insKey, lupKey;

int main(int argc, char **argv)
{
	// Parse the command line
	char c;
	int seed  = 0;
	int cmd_n = 100000;
	int cmd_q = 1000000;
    int times = 1;
	FILE *fp  = stdin;
	while ((c = getopt(argc, argv, "r:f:n:q:t:")) != EOF) {
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
            
		case 't':
			times = atoi(optarg);
            printf("times == %d\n", times);
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

    test_all();
    //test_memory_usage();

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
    vector<double> t1, t2;
    printf("maxSteps %d\n", maxSteps);
    double sum1 = 0;
    double sum2 = 0;

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

        StandardCuckooFilter<uint16_t, 16> sscf;
        m = 1 << (int)(ceil(log2(n / b)));
        sscf.init(m, b, maxSteps);
        double sb = evaluate(sscf, "Semi-Sort-Standard-CuckooFilter", insKey, lupKey);
        t1.push_back(sb);
        printf("%.5f\n", sb);
        sum1 += sb;

        VacuumFilter<uint16_t, 16> xor_filter;
        //m = 1 << (int)(ceil(log2(n / b / 0.95)));
        m = (int)(n / b);
        //m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        xor_filter.init(m, b, maxSteps);
        //xor_filter.test_bucket();
        sb = evaluate(xor_filter, "Vacuum-CuckooFilter", insKey, lupKey);
        t2.push_back(sb);
        printf("%.5f\n", sb);
        sum2 += sb;

        RandomFilter<uint16_t, 8> xor_filter;
        //m = 1 << (int)(ceil(log2(n / b / 0.95)));
        m = n / b;
        //m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        xor_filter.init(m, b, maxSteps);
        double sb = evaluate(xor_filter, "random-CuckooFilter", insKey, lupKey);
        printf("%.5f\n", sb);
        sum1 += sb;

        AddSubFilter<uint8_t, 8> new_filter;
        //m = 1 << (int)(ceil(log2(n / b / 0.95)));
        //m = (int) (n / b / 0.96 + 1);
        m = n / b;
        m += m & 1;
        new_filter.init(m, b, maxSteps);
        sum2 += evaluate(new_filter, "addsub-CuckooFilter", insKey, lupKey);
    }

    printf("%.5f %.5f\n", sum1 / times, sum2 / times);
    */

    /*
    sort(t1.begin(), t1.end());
    sort(t2.begin(), t2.end());

    printf("t1 : %.5f %.5f\n", t1[0], t1[1]);
    printf("t2 : %.5f %.5f\n", t2[0], t2[1]);

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

    /*
    MortonAddFilter<uint8_t, 8> new_filter;
    //m = 1 << (int)(ceil(log2(n / b / 0.95)));
    m = (int) (n / b / 0.95 + 1);
    //m = n / b;
    m += m & 1;
    new_filter.init(m, b, maxSteps);
    evaluate(new_filter, "Morton-CuckooFilter", insKey, lupKey);

    VacuumFilter<uint16_t, 12> vac_filter;
    //m = 1 << (int)(ceil(log2(n / b / 0.95)));
    //m = n / b;
    m = (int) (n / b / 0.96 + 1);
    m += m & 1;
    vac_filter.init(m, b, maxSteps);
    evaluate(vac_filter, "vacuume-CuckooFilter", insKey, lupKey);
    //sum1 += sb;
    */

    /*
    BloomFilter<uint8_t, 8> bloom_filter;
	m = (int) (n / b / 0.96 + 1);
	m += m & 1;
	bloom_filter.init(m, b, maxSteps);
	evaluate(bloom_filter, "bloom-CuckooFilter", insKey, lupKey);

    AddSubFilter<uint16_t, 8> add_filter;
    //m = 1 << (int)(ceil(log2(n / b / 0.95)));
    //m = n / b;
    m = (int) (n / b / 0.96 + 1);
    m += m & 1;
    add_filter.init(m, b, maxSteps);
    evaluate(add_filter, "add-sub-CuckooFilter", insKey, lupKey);
    printf("full factor = %.5f\n", add_filter.get_full_bucket_factor());
    //sum1 += sb;
    */

	return 0;	
}
