#include <bits/stdc++.h>
#include <random>
#include <unistd.h>
#include "hashutil.h"
#include "cuckoo.h"
#include <time.h>
#include <ratio>
#include <chrono>
#include <cuckoofilter/cuckoofilter.h>

#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a)  cerr << #a << " = " << a << endl
#define memcle(a) memset(a, 0, sizeof(a))

using namespace std;

template <typename fp_t, int fin_len>
void evaluate(Filter<fp_t, fin_len> &filter,
			  const char *filter_name,
			  const vector<int> &insKey,
			  const vector<int> &lupKey,
              vector<double> &result,
              bool verbose = false)
{
	int fail_insert = 0;
	int false_positive = 0;
	int false_negative = 0;

    printf("insert\n");
	
	for (int key : insKey)
    {
		if (filter.insert(key) == 1) // insertion failed
        {
            //printf("insert false : %d\n", key);
			++fail_insert;
            //filter.insert(key);
            break;
        }
        //printf("%.5f\n", filter.get_load_factor());
    }

    /*
    printf("positive lookup\n");
	for (int key : insKey)
    {
		if (filter.lookup(key) == false) // false negative
        {
            //filter.lookup(key);
            //printf("lookup false : %d\n", key);
            //filter.lookup(key);
			++false_negative;
        }
    }
    */

    printf("negative lookup\n");

	//unordered_set<int> S(insKey.begin(), insKey.end());

    //printf("cannot generate a set");

	int n = insKey.size(), q = lupKey.size(), m = filter.n, b = filter.m;

    deln(q);
	for (int i = 0; i < q; i++)
    {
        int key = lupKey[i];
		//if (filter.lookup(key) == true && !S.count(key)) // false positive
        if (filter.lookup(key))
			++false_positive;
    }

    printf("negative lookup finish\n");


    if (verbose)
    {
        printf("testing %s: n = %d, q = %d, m = %d, b = %d memory = %lu bytes\n",
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
    out = fopen("key-memory-low-fp.csv", "w");
    assert(out != NULL);

    fprintf(out, "key number, Semi-Sort-Cuckoo-Filter-16bit, Morton-Filter, Vacuum-Filter-16bit, Vacuum-Filter-16bit, Bloom-Filter-20bit\n");
    int max = 10000000;

    for (int n = 1000; n <= max; n += 1000)
    {
        if (n % (max / 100) == 0)
            printf("process : %d%%\n", n / (max / 100));

        StandardCuckooFilter<uint16_t, 16> sscf;
        m = 1 << (int)(ceil(log2(n / b / 0.95)));
        sscf.init(m, b, maxSteps);

        MortonAddFilter<uint8_t, 8> mf;
        m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        mf.init(m, b, maxSteps);

        VacuumFilter<uint16_t, 16> vf;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        vf.init(m, b, maxSteps);

        VacuumFilter<uint16_t, 16> vf2;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        vf2.init(m, b, maxSteps);

        BloomFilter<uint16_t, 20> bf;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        bf.init(m, b, maxSteps);

        fprintf(out, "%d, %lu, %lu, %lu, %lu, %lu\n", n, sscf.memory_consumption, mf.memory_consumption, vf.memory_consumption, vf2.memory_consumption, bf.memory_consumption);
    }

    fclose(out);
}

void test_all()
{
    FILE *out = NULL;
    out = fopen("all-performance.csv", "w");
    assert(out != NULL);

    //result : 
    //0 : memory usage
    //1 : fail insert rate
    //2 : false negative rate
    //3 : false positive rate
    //4 : load factor
    //5 : bits per item

    VacuumFilter<uint16_t, 13> sb;
    sb.init(100, 4, 200);
    //sb.test_bucket();

    printf("pass\n");
    fprintf(out, "key size, ");

    //const int test_time = 10;
    //const int key[test_time] = {100000, 200000, 300000, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, int((1 << 20) * 0.94)};
    const int test_time = 3;
    const int key[test_time] = {1000000, 15000000, 20000000};
    const string name[5] = {"Vacuum", "Standard", "Morton", "Bloom", "Offset"};
    const string metric[6] = {"memory", "fail insert", "false negative", "false positive", "load factor", "bits per item"};

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 5; j++)
        {
            string tmp = name[j] + " " + metric[i];
            if (i != 5 || j != 4)
                fprintf(out, "%s, ", tmp.c_str());
            else 
                fprintf(out, "%s\n", tmp.c_str());
        }

    for (int i = 0; i < 1; i++)
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
            //insKey.push_back(rd());
            insKey.push_back(i);
        for (int i = 1; i <= q; ++i)
            lupKey.push_back(i + n);

        vector<double> res[5];

        MortonAddFilter<uint8_t, 8> mf;
        m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        mf.init(m, b, maxSteps);
        evaluate(mf, name[2].c_str(),  insKey, lupKey, res[2]);
        printf("mf finish\n");

        VacuumFilter<uint16_t, 13> vf;
        m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        vf.init(m, b, maxSteps);
        evaluate(vf, name[0].c_str(), insKey, lupKey, res[0]);

        printf("vf finish\n");

        StandardCuckooFilter<uint16_t, 13> cf;
        m = 1 << (int)(ceil(log2(n / b / 0.95)));
        cf.init(m, b, maxSteps);
        evaluate(cf, name[1].c_str(),  insKey, lupKey, res[1]);

        printf("cf finish\n");



        BloomFilter<uint16_t, 14> bf;
        m = (int) (n / b / 0.95 + 1);
        m += m & 1;
        bf.init(m, b, maxSteps);
        evaluate(bf, name[3].c_str(),  insKey, lupKey, res[3]);

        printf("bf finish\n");

        OffsetFilter<uint16_t, 12> of;
        m = (int) (n / b / 0.96 + 1);
        m += m & 1;
        of.init(m, b, maxSteps);
        evaluate(of, name[4].c_str(), insKey, lupKey, res[4]);

        printf("of finish\n");

        fprintf(out, "%d, ", n);

        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 5; k++)
            {
                if (j != 5 || k != 4)
                    fprintf(out, "%.5f, ", res[k][j]);
                else
                    fprintf(out, "%.5f\n", res[k][j]);
            }
    }
}

void test_load_factor_new()
{
    FILE *out = NULL;
    out = fopen("load-factor-new.csv", "w");
    assert(out != NULL);

    const int test_time = 7;
    static constexpr int alt_len[7] = {64, 128, 256, 512, 1024, 2048, 4096};

    int n = 1 << 27;
    int q = 1;
    int maxSteps = 400;
    int b = 4;
    int seed = 1;
    int m;
    int rept = 1;
    const int fp_len = 16;
    mt19937 rd(seed);


    vector<int> insKey;
    vector<int> lupKey;

    for (int i = 1; i <= n; ++i)
        insKey.push_back(rd());
    for (int i = 1; i <= q; ++i)
        lupKey.push_back(rd());

    cuckoofilter::VacuumFilter<size_t, 12> vf(n * 0.95);
    for (int i = 0; i < n; i++)
        if (vf.Add(insKey[i]) != cuckoofilter::Ok)
        {
            printf("%.5f\n", vf.LoadFactor());
            break;
        }

    return;

    fprintf(out, "alt_len, vf load factor, item numbers = %d, max kick = %d\n", n, maxSteps);
    for (int i = 0; i < test_time; i++)
    {
        double lf = 0;
        for (int j = 0; j < rept; j++)
        {
            vector<int> insKey;
            vector<int> lupKey;

            for (int i = 1; i <= n; ++i)
                insKey.push_back(rd());
            for (int i = 1; i <= q; ++i)
                lupKey.push_back(rd());

            vector<double> res[4];

            VacuumFilterTest<uint16_t, fp_len> vf;
            m = int(n / b);
            m += m & 1;
            vf.init(m, b, maxSteps);
            vf.set_alt_len(alt_len[i]);
            evaluate(vf, "vf", insKey, lupKey, res[0]);

            lf += res[0][4] / rept;
        }

        printf("%d, %.5f\n", alt_len[i], lf);
        fprintf(out, "%d, %.5f\n", alt_len[i], lf);
    }
}

void test_load_factor()
{
    FILE *out = NULL;
    out = fopen("load-factor.csv", "w");
    assert(out != NULL);

    //result : 
    //0 : memory usage
    //1 : fail insert rate
    //2 : false negative rate
    //3 : false positive rate
    //4 : load factor
    //5 : bits per item

    fprintf(out, "key size, ");

    const int test_time = 14;
    //const int key[test_time] = {5000000, 6000000, 7000000, 8000000};
    const int key[test_time] = {int((1 << 27) * 0.95), 100000, 200000, 300000, 400000, 500000, 1000000, 2000000, 10000000, 4000000,5000000, 6000000, 7000000, 8000000};
    const string name[4] = {"Vacuum", "Standard", "Morton", "Offset"};
    const string metric[6] = {"memory", "fail insert", "false negative", "false positive", "load factor", "bits per item"};

    for (int j = 0; j < 4; j++)
        {
            int i = 4;
            string tmp = name[j] + " " + metric[i];
            if (j != 4 - 1)
                fprintf(out, "%s, ", tmp.c_str());
            else
                fprintf(out, "%s\n", tmp.c_str());
        }

    for (int i = 0; i < 1; i++)
    {
        int n = key[i];
        int q = 1;
        int maxSteps = 400;
        int b = 4;
        int seed = 1;
        int m;
        int rept = 5;
        mt19937 rd(seed);

        double load_factor[4];
        memset(load_factor, 0, sizeof(load_factor));

        for (int j = 0; j < rept; j++)
        {
            vector<int> insKey;
            vector<int> lupKey;

            for (int i = 1; i <= n; ++i)
                insKey.push_back(rd());
            for (int i = 1; i <= q; ++i)
                lupKey.push_back(rd());

            vector<double> res[4];

            VacuumFilter<uint16_t, 16> vf;
            m = int(n / b);
            m += m & 1;
            vf.init(m, b, maxSteps);
            evaluate(vf, name[0].c_str(), insKey, lupKey, res[0]);

            /*
            StandardCuckooFilter<uint16_t, 16> cf;
            m = 1 << (int)(ceil(log2(n / b)));
            cf.init(m, b, maxSteps);
            evaluate(cf, name[1].c_str(),  insKey, lupKey, res[1]);

            MortonAddFilter<uint8_t, 8> mf;
            m = (int) (n / b);
            m += m & 1;
            mf.init(m, b, maxSteps);
            evaluate(mf, name[2].c_str(),  insKey, lupKey, res[2]);

            OffsetFilter<uint16_t, 16> of;
            m = int(n / b);
            m += m & 1;
            of.init(m, b, maxSteps);
            evaluate(of, name[4].c_str(), insKey, lupKey, res[3]);
            */

            for (int k = 0; k < 1; k++)
                load_factor[k] += res[k][4] / rept;
        }

        fprintf(out, "%d, ", n);

        for (int k = 0; k < 4; k++)
        {
            if (k != 3)
                fprintf(out, "%.5f, ", load_factor[k]);
            else
                fprintf(out, "%.5f\n", load_factor[k]);
        }
    }
}

double time_cost(chrono::steady_clock::time_point &start,
                 chrono::steady_clock::time_point &end)
{
     /* Return the time elapse between start and end
      * count by seconds
      */
    double elapsedSeconds = ((end - start).count()) * chrono::steady_clock::period::num / static_cast<double>(chrono::steady_clock::period::den);
    return elapsedSeconds;
}

void test_ins()
{
    FILE *out = fopen("throughput_ins.csv", "w");
    assert(out != NULL);

    int rept = 1;
    int n = (1 << 20) * 0.94;
    int q = 1000000;
    int maxSteps = 400;
    int b = 4;
    int seed = 1;
    int m;
    static constexpr int fin_len[6] = {13, 14, 13, 13, 12, 12};

    mt19937 rd(seed);
    double mop[6][30], mop1[6][30];
    int cnt[6][30], cnt1[6][30];

    memcle(mop);
    memcle(cnt);
    memcle(mop1);
    memcle(cnt1);

    printf("vacuum insert\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        VacuumFilter<uint16_t, fin_len[0]> vf;
        m = n / b;
        vf.init(m, b, maxSteps);

        int i = 0;
        int j = 0;
        int lim;

        auto tot_start = chrono::steady_clock::now();

        for (double r = 0.05; r <= 0.96; r += 0.05, ++j)
        {
            //printf("%.2f\n", r);
            lim = int(n * r);
            int ins_cnt = 0;

            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                //printf("%d\n", i);
                vf.insert(insKey[i]);
                ++ins_cnt;
            }

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[0][j] += double(ins_cnt) / 1000000.0 / cost;
            cnt[0][j] += 1;
            //fprintf(out, "%d, %.2f, %.5f, %.5f\n", lim, r, cost, double(ins_cnt) / 1000000.0 / cost);
        }

        auto tot_end = chrono::steady_clock::now();

        mop[0][j] += double(lim) / 1000000.0 / time_cost(tot_start, tot_end);
        cnt[0][j] += 1;
    }


    printf("bf insert\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        BloomFilter<uint16_t, fin_len[1]> vf;
        m = n / b;
        vf.init(m, b, maxSteps);

        int i = 0;
        int j = 0;
        int lim;

        auto tot_start = chrono::steady_clock::now();

        for (double r = 0.05; r <= 0.96; r += 0.05, ++j)
        {
            //printf("%.2f\n", r);
            lim = int(n * r);
            int ins_cnt = 0;

            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                //printf("%d\n", i);
                vf.insert(insKey[i]);
                ++ins_cnt;
            }

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[1][j] += double(ins_cnt) / 1000000.0 / cost;
            cnt[1][j] += 1;
            //fprintf(out, "%d, %.2f, %.5f, %.5f\n", lim, r, cost, double(ins_cnt) / 1000000.0 / cost);
        }

        auto tot_end = chrono::steady_clock::now();

        mop[1][j] += double(lim) / 1000000.0 / time_cost(tot_start, tot_end);
        cnt[1][j] += 1;
    }

    printf("fb-vacuum-ss insert\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::VacuumFilter<size_t, fin_len[2], cuckoofilter::PackedTable > cf(n);
        int i = 0;
        int j = 0;
        int lim;

        auto tot_start = chrono::steady_clock::now();

        for (double r = 0.05; r <= 0.96; r += 0.05, ++j)
        {
            //printf("%.2f\n", r);
            lim = int(n * r);
            int ins_cnt = 0;

            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                cf.Add(insKey[i]);
                ++ins_cnt;
            }

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[2][j] += double(ins_cnt) / 1000000.0 / cost;
            cnt[2][j] += 1;
            //fprintf(out, "%d, %.2f, %.5f, %.5f\n", lim, r, cost, double(ins_cnt) / 1000000.0 / cost);
        }

        auto tot_end = chrono::steady_clock::now();

        mop[2][j] += double(lim) / 1000000.0 / time_cost(tot_start, tot_end);
        cnt[2][j] += 1;
     }

    printf("fb-std-ss insert\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::CuckooFilter<size_t, fin_len[3], cuckoofilter::PackedTable > cf(n);
        int i = 0;
        int j = 0;
        int lim;

        auto tot_start = chrono::steady_clock::now();

        for (double r = 0.05; r <= 0.96; r += 0.05, ++j)
        {
            //printf("%.2f\n", r);
            lim = int(n * r);
            int ins_cnt = 0;

            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                cf.Add(insKey[i]);
                ++ins_cnt;
            }

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[3][j] += double(ins_cnt) / 1000000.0 / cost;
            cnt[3][j] += 1;
            //fprintf(out, "%d, %.2f, %.5f, %.5f\n", lim, r, cost, double(ins_cnt) / 1000000.0 / cost);
        }

        auto tot_end = chrono::steady_clock::now();

        mop[3][j] += double(lim) / 1000000.0 / time_cost(tot_start, tot_end);
        cnt[3][j] += 1;
     }

    printf("fb-vacuum insert\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::VacuumFilter<size_t, fin_len[4]> cf(n);
        int i = 0;
        int j = 0;
        int lim;

        auto tot_start = chrono::steady_clock::now();

        for (double r = 0.05; r <= 0.96; r += 0.05, ++j)
        {
            //printf("%.2f\n", r);
            lim = int(n * r);
            int ins_cnt = 0;

            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                cf.Add(insKey[i]);
                ++ins_cnt;
            }

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[4][j] += double(ins_cnt) / 1000000.0 / cost;
            cnt[4][j] += 1;
            //fprintf(out, "%d, %.2f, %.5f, %.5f\n", lim, r, cost, double(ins_cnt) / 1000000.0 / cost);
        }

        auto tot_end = chrono::steady_clock::now();

        mop[4][j] += double(lim) / 1000000.0 / time_cost(tot_start, tot_end);
        cnt[4][j] += 1;
    }

    printf("fb-std insert\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::CuckooFilter<size_t, fin_len[5]> cf(n);
       int i = 0;
        int j = 0;
        int lim;

        auto tot_start = chrono::steady_clock::now();

        for (double r = 0.05; r <= 0.96; r += 0.05, ++j)
        {
            //printf("%.2f\n", r);
            lim = int(n * r);
            int ins_cnt = 0;

            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                cf.Add(insKey[i]);
                ++ins_cnt;
            }

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[5][j] += double(ins_cnt) / 1000000.0 / cost;
            cnt[5][j] += 1;
            //fprintf(out, "%d, %.2f, %.5f, %.5f\n", lim, r, cost, double(ins_cnt) / 1000000.0 / cost);
        }

        auto tot_end = chrono::steady_clock::now();

        mop[5][j] += double(lim) / 1000000.0 / time_cost(tot_start, tot_end);
        cnt[5][j] += 1;
    }


    fprintf(out, "occupancy, VF-ss-no-%d, BF-%d, VF-ss-%d, CF-ss-%d, VF-%d, CF-%d(MOPS), key-size-%d\n",  fin_len[0], fin_len[1], fin_len[2], fin_len[3], fin_len[4], fin_len[5], n);

    for (int j = 0; j < 20; j++) 
    {
        fprintf(out, "%.2f, ", (j + 1) * 0.05);
        for (int k = 0; k < 6; k++)
            fprintf(out, "%.5f, ", mop[k][j] / cnt[k][j]);
        fprintf(out, "\n");
    }
}


void test_lookup()
{

    FILE *out = fopen("throughput_lookup.csv", "w");
    assert(out != NULL);

    int rept = 1;
    int n = (1 << 17) * 0.94;
    int q = 100000;
    int maxSteps = 400;
    int b = 4;
    int seed = 1;
    int m;
    double neg_frac = 0.5;
    static constexpr int fin_len[6] = {13, 14, 13, 13, 12, 12};

    mt19937 rd(seed);
    double mop[6][20], mop1[6][20];
    int cnt[6][20], cnt1[6][20];
    memcle(mop);
    memcle(cnt);
    memcle(mop1);
    memcle(cnt1);

    printf("vf lookup\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        VacuumFilter<uint16_t, fin_len[0]> vf;
        m = n / b;
        vf.init(m, b, maxSteps);

        int i = 0;
        int ins_cnt = 0;
        int j = 0;

        for (double r = 0.05; r <= 0.96; r += 0.05, j++)
        {
            printf("%.2f\n", r);
            int lim = int(n * r);
            for (; i < lim; i++) 
                if (vf.insert(insKey[i]) == 0)
                    ++ins_cnt;

            auto start = chrono::steady_clock::now();

            int lookup_number = 0;
            int p = min(q, lim);
            int k = 0;

            for (k = 0; k < p * neg_frac; k++)
                if (vf.lookup(lupKey[k]) == 0)
                    lookup_number++;

            for (; k < p; k++)
                if (vf.lookup(insKey[k]) == 1)
                    lookup_number++;

            auto end = chrono::steady_clock::now();

            double cost = time_cost(start, end);
            mop[0][j] += double(p) / 1000000.0 / cost;
            cnt[0][j] += 1;

            start = chrono::steady_clock::now();
            int t = min(q, lim);
            for (int i = 0; i < t; i++)
                vf.lookup(insKey[i]);

            end = chrono::steady_clock::now();

            cost = time_cost(start, end);
            mop1[0][j] += double(t) / 1000000.0 / cost;
            cnt1[0][j] += 1;
            //fprintf(out, "%.5f, %.5f\n", cost, double(t) / 1000000.0 / cost);
        }
    }

    printf("bloom lookup\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        BloomFilter<uint16_t, fin_len[1]> vf;
        m = n / b;
        vf.init(m, b, maxSteps);

        int i = 0;
        int ins_cnt = 0;
        int j = 0;

        for (double r = 0.05; r <= 0.96; r += 0.05, j++)
        {
            printf("%.2f\n", r);
            int lim = int(n * r);
            for (; i < lim; i++) 
                if (vf.insert(insKey[i]) == 0)
                    ++ins_cnt;

            auto start = chrono::steady_clock::now();

            int lookup_number = 0;
            int p = min(q, lim);
            int k = 0;

            for (k = 0; k < p * neg_frac; k++)
                if (vf.lookup(lupKey[k]) == 0)
                    lookup_number++;

            for (; k < p; k++)
                if (vf.lookup(insKey[k]) == 1)
                    lookup_number++;

            auto end = chrono::steady_clock::now();

            double cost = time_cost(start, end);
            mop[1][j] += double(p) / 1000000.0 / cost;
            cnt[1][j] += 1;

            start = chrono::steady_clock::now();
            int t = min(q, lim);
            for (int i = 0; i < t; i++)
                vf.lookup(insKey[i]);

            end = chrono::steady_clock::now();

            cost = time_cost(start, end);
            mop1[1][j] += double(t) / 1000000.0 / cost;
            cnt1[1][j] += 1;
            //fprintf(out, "%.5f, %.5f\n", cost, double(t) / 1000000.0 / cost);
        }
    }

    printf("fb-vacuum-ss lookup\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());
       
        cuckoofilter::VacuumFilter<size_t, fin_len[2], cuckoofilter::PackedTable> vf(n);

        int i = 0;
        int ins_cnt = 0;
        int j = 0;

        for (double r = 0.05; r <= 0.96; r += 0.05, j++)
        {
            printf("%.2f\n", r);
            int lim = int(n * r);
            for (; i < lim; i++) 
                if (vf.Add(insKey[i]) == 0)
                    ++ins_cnt;

            auto start = chrono::steady_clock::now();

            int lookup_number = 0;
            int p = min(q, lim);
            int k = 0;

            for (k = 0; k < p * neg_frac; k++)
                if (vf.Contain(lupKey[k]) != cuckoofilter::Ok)
                    lookup_number++;

            for (; k < p; k++)
                if (vf.Contain(insKey[k]) == cuckoofilter::Ok)
                    lookup_number++;

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);

            mop[2][j] += double(p) / 1000000.0 / cost;
            cnt[2][j] += 1;


            start = chrono::steady_clock::now();

            int t = min(q, lim);
            for (int i = 0; i < t; i++)
                if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                    lookup_number++;

            end = chrono::steady_clock::now();
            cost = time_cost(start, end);
            mop1[2][j] += double(t) / 1000000.0 / cost;
            cnt1[2][j] += 1;

            printf("%d\n", lookup_number);
        }
    }

    printf("fb-std-ss lookup\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::CuckooFilter<size_t, fin_len[3], cuckoofilter::PackedTable> vf(n);
        int i = 0;
        int ins_cnt = 0;
        int j = 0;

        for (double r = 0.05; r <= 0.96; r += 0.05, j++)
        {
            printf("%.2f\n", r);
            int lim = int(n * r);
            for (; i < lim; i++) 
                if (vf.Add(insKey[i]) == 0)
                    ++ins_cnt;

            auto start = chrono::steady_clock::now();

            int lookup_number = 0;
            int p = min(q, lim);
            int k = 0;

            for (k = 0; k < p * neg_frac; k++)
                if (vf.Contain(lupKey[k]) != cuckoofilter::Ok)
                    lookup_number++;

            for (; k < p; k++)
                if (vf.Contain(insKey[k]) == cuckoofilter::Ok)
                    lookup_number++;

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);

            mop[3][j] += double(p) / 1000000.0 / cost;
            cnt[3][j] += 1;

            start = chrono::steady_clock::now();

            int t = min(q, lim);
            for (int i = 0; i < t; i++)
                if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                    lookup_number++;

            end = chrono::steady_clock::now();
            cost = time_cost(start, end);
            mop1[3][j] += double(t) / 1000000.0 / cost;
            cnt1[3][j] += 1;

            printf("%d\n", lookup_number);
        }
    }

    printf("fb-vacuum lookup\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::VacuumFilter<size_t, fin_len[4]> vf(n);

        int i = 0;
        int ins_cnt = 0;
        int j = 0;

        for (double r = 0.05; r <= 0.96; r += 0.05, j++)
        {
            printf("%.2f\n", r);
            int lim = int(n * r);
            for (; i < lim; i++) 
                if (vf.Add(insKey[i]) == 0)
                    ++ins_cnt;

            auto start = chrono::steady_clock::now();

            int lookup_number = 0;
            int p = min(q, lim);
            int k = 0;

            for (k = 0; k < p * neg_frac; k++)
                if (vf.Contain(lupKey[k]) != cuckoofilter::Ok)
                    lookup_number++;

            for (; k < p; k++)
                if (vf.Contain(insKey[k]) == cuckoofilter::Ok)
                    lookup_number++;

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);

            mop[4][j] += double(p) / 1000000.0 / cost;
            cnt[4][j] += 1;

            start = chrono::steady_clock::now();

            int t = min(q, lim);
            for (int i = 0; i < t; i++)
                if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                    lookup_number++;

            end = chrono::steady_clock::now();
            cost = time_cost(start, end);
            mop1[4][j] += double(t) / 1000000.0 / cost;
            cnt1[4][j] += 1;

            printf("%d\n", lookup_number);
        }
    }

    printf("fb-std lookup\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::CuckooFilter<size_t, fin_len[5]> vf(n);
        
        int i = 0;
        int ins_cnt = 0;
        int j = 0;

        for (double r = 0.05; r <= 0.96; r += 0.05, j++)
        {
            printf("%.2f\n", r);
            int lim = int(n * r);
            for (; i < lim; i++) 
                if (vf.Add(insKey[i]) == 0)
                    ++ins_cnt;

            auto start = chrono::steady_clock::now();

            int lookup_number = 0;
            int p = min(q, lim);
            int k = 0;

            for (k = 0; k < p * neg_frac; k++)
                if (vf.Contain(lupKey[k]) != cuckoofilter::Ok)
                    lookup_number++;

            for (; k < p; k++)
                if (vf.Contain(insKey[k]) == cuckoofilter::Ok)
                    lookup_number++;

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);

            mop[5][j] += double(p) / 1000000.0 / cost;
            cnt[5][j] += 1;

            start = chrono::steady_clock::now();

            int t = min(q, lim);
            for (int i = 0; i < t; i++)
                if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                    lookup_number++;

            end = chrono::steady_clock::now();
            cost = time_cost(start, end);
            mop1[5][j] += double(t) / 1000000.0 / cost;
            cnt1[5][j] += 1;

            printf("%d\n", lookup_number);
        }
     }

    fprintf(out, "occupancy, VF-ss-no-%d neg, VF-ss-no pos, BF-%d neg, BF pos, VF-ss-%d neg, VF-ss pos, CF-ss-%d neg, CF-ss pos, VF-%d neg, VF pos, CF-%d neg, CF pos, negative fraction = %.2f, item numbers = %d, query number = %d\n", fin_len[0], fin_len[1], fin_len[2], fin_len[3], fin_len[4], fin_len[5], neg_frac, n, q);

    for (int j = 0; j < 19; j++)
    {
        fprintf(out, "%.2f, ", (j + 1) * 0.05);
        for (int k = 0; k < 6; k++)
            fprintf(out, "%.5f, %.5f, ", mop[k][j] / cnt[k][j], mop1[k][j] / cnt1[k][j]);
        fprintf(out, "\n");
    }
}

void test_del()
{
    FILE *out = fopen("throughput_del.csv", "w");
    assert(out != NULL);

    int rept = 1;
    int n = (1 << 17) * 0.94;
    int q = 1;
    int maxSteps = 400;
    int b = 4;
    int seed = 1;
    int m;
    static constexpr int fin_len[5] = {13, 13, 13, 12, 12};

    mt19937 rd(seed);
    double mop[6][20], mop1[6][20];
    int cnt[6][20], cnt1[6][20];
    memcle(mop);
    memcle(cnt);
    memcle(mop1);
    memcle(cnt1);
    
    printf("vf-ss-no delete\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        VacuumFilter<uint16_t, fin_len[0]> vf;
        m = n / b;
        vf.init(m, b, maxSteps);

        for (int i = 0; i < 0.95 * n; i++) vf.insert(insKey[i]);

        int i = 0;
        int j = 18;
        for (double r = 0.05; r <= 0.96; r += 0.05, --j)
        {
            int del_cnt = 0;
            int lim = int(n * r);
            //deln(lim);
            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                if (vf.del(insKey[i]) == true)
                    ++del_cnt;
            }

            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);

            mop[0][j] += double(del_cnt) / 1000000.0 / cost;
            cnt[0][j] += 1;
        }

        /*
        for (int i = 0; i < 0.95 * n; i++)
            if (vf.lookup(insKey[i]))
                printf("i = %d\n", i);
        */
    }

    printf("vf-ss delete\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::VacuumFilter<size_t, fin_len[1], cuckoofilter::PackedTable> vf(n);

        for (int i = 0; i < 0.95 * n; i++) vf.Add(insKey[i]);

        int i = 0;
        int j = 18;
        for (double r = 0.05; r <= 0.96; r += 0.05, --j)
        {
            int del_cnt = 0;
            int lim = int(n * r);
            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                if (vf.Delete(insKey[i]) == cuckoofilter::Ok)
                    ++del_cnt;
            }
            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[1][j] += double(del_cnt) / 1000000.0 / cost;
            cnt[1][j] += 1;
        }

        /*
        for (int i = 0; i < 0.95 * n; i++)
            if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                printf("i = %d\n", i);
        */
    }

    printf("cf-ss delete\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::CuckooFilter<size_t, fin_len[2], cuckoofilter::PackedTable> vf(n);

        for (int i = 0; i < 0.95 * n; i++) vf.Add(insKey[i]);

        int i = 0;
        int j = 18;
        for (double r = 0.05; r <= 0.96; r += 0.05, --j)
        {
            int del_cnt = 0;
            int lim = int(n * r);
            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                if (vf.Delete(insKey[i]) == cuckoofilter::Ok)
                    ++del_cnt;
            }
            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[2][j] += double(del_cnt) / 1000000.0 / cost;
            cnt[2][j] += 1;
        }

        /*
        for (int i = 0; i < 0.95 * n; i++)
            if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                printf("i = %d\n", i);
        */
    }

    printf("vf delete\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::VacuumFilter<size_t, fin_len[3]> vf(n);

        for (int i = 0; i < 0.95 * n; i++) vf.Add(insKey[i]);

        int i = 0;
        int j = 18;
        for (double r = 0.05; r <= 0.96; r += 0.05, --j)
        {
            int del_cnt = 0;
            int lim = int(n * r);
            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                if (vf.Delete(insKey[i]) == cuckoofilter::Ok)
                    ++del_cnt;
            }
            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);
            mop[3][j] += double(del_cnt) / 1000000.0 / cost;
            printf("%.2f\n", double(del_cnt) / 1000000.0 / cost);
            cnt[3][j] += 1;
        }

        /*
        for (int i = 0; i < 0.95 * n; i++)
            if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                printf("i = %d\n", i);
        */
    }

    printf("cf delete\n");
    for (int t = 0; t < rept; t++)
    {
        vector<int> insKey;
        vector<int> lupKey;
        for (int i = 1; i <= n; ++i) insKey.push_back(rd());
        for (int i = 1; i <= q; ++i) lupKey.push_back(rd());

        cuckoofilter::CuckooFilter<size_t, 12> vf(n);

        for (int i = 0; i < 0.95 * n; i++) vf.Add(insKey[i]);

        int i = 0;
        int j = 18;
        for (double r = 0.05; r <= 0.96; r += 0.05, --j)
        {
            int del_cnt = 0;
            int lim = int(n * r);
            auto start = chrono::steady_clock::now();
            for (; i < lim; i++) 
            {
                if (vf.Delete(insKey[i]) == cuckoofilter::Ok)
                    ++del_cnt;
            }
            auto end = chrono::steady_clock::now();
            double cost = time_cost(start, end);

            /*
            printf("time = %.2f\n", cost);
            printf("del = %d\n", del_cnt);
            printf("mops = %.2f\n", double(del_cnt) / 1000000.0 / cost);
            */

            printf("%.2f\n", double(del_cnt) / 1000000.0 / cost);
            mop[4][j] += double(del_cnt) / 1000000.0 / cost;
            cnt[4][j] += 1;
        }

        /*
        for (int i = 0; i < 0.95 * n; i++)
            if (vf.Contain(insKey[i]) == cuckoofilter::Ok)
                printf("i = %d\n", i);
        */
    }

    fprintf(out, "occupancy, VF-ss-no-%d del, VF-ss-%d del, CF-ss-%d del, VF-%d del, CF-%d del, item numbers = %d\n", fin_len[0], fin_len[1], fin_len[2], fin_len[3], fin_len[4], n);
    for (int j = 0; j < 19; j++) 
    {
        fprintf(out, "%.2f, ", (j + 1) * 0.05);
        for (int k = 0; k < 5; k++)
            fprintf(out, "%.5f, ", mop[k][j] / cnt[k][j]);
        fprintf(out, "\n");
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
	int cmd_q = 8000000;
    string test_name;
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
			//times = atoi(optarg);
            test_name = optarg;
            printf("test_name == %s\n", test_name.c_str());
			break;
            
		default:
			break;
		}	
	}

	int n, q; // #item, #lookup, #bucket, #entry per bucket

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
        deln(q);
		for (int i = 1; i <= n; ++i)
			insKey.push_back(rd());
		for (int i = 1; i <= q; ++i)
			lupKey.push_back(rd());
	}
    
    printf("testing with parameters : \n");
    printf("insert keys : %d\n", n);
    printf("lookup times: %d\n\n\n", q);

    if (test_name == "all") test_all();
    if (test_name == "memory") test_memory_usage();
    if (test_name == "load") test_load_factor();
    if (test_name == "ins") test_ins();
    if (test_name == "lookup") test_lookup();
    if (test_name == "del") test_del();
    if (test_name == "newvf") test_load_factor_new();

    /*
	MortonAddFilter<uint16_t, 8> mf;
	//m = 1 << (int)(ceil(log2(n / b / 0.95)));
	m = (int) (n / b / 0.95 + 1);
	m += m & 1;
	mf.init(m, b, maxSteps);
	evaluate(mf, "MortonFilter", insKey, lupKey);
    */

    /*
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
