#include <bits/stdc++.h>
#include <iostream>
#include <chrono>
#include <ctime>
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#define deln(x) cerr << #x << " " << x << endl
using namespace std;

inline int find_highest_bit(int v)
{
// tricks of bit 
// from http://graphics.stanford.edu/~seander/bithacks.html
    static const int MultiplyDeBruijnBitPosition[32] = 
    {
      0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
      8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };

    v |= v >> 1; // first round down to one less than a power of 2 
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    int r = MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
    return r;
}


unsigned long lower_power_of_two(unsigned long x)
{
    if (x >= 64) return 64;
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    //x = x | (x >> 8);
    //x = x | (x >> 16);
    return x - (x >> 1);
}

int n = 1 << 27;
int max_2_power = 1 << 27;

int benchmark(int pos, int fp) // get alternate position
{

    int fp_hash = fp * 0x5bd1e995;

    int delta = fp_hash & (max_2_power - 1);

    return pos ^ delta;

    int bias = delta;

    pos = pos + bias;
    if (pos >= n) pos -= n;
    
    int segment_length = 1 << find_highest_bit(pos ^ n);
    int t = (delta & (segment_length - 1) & 255);
    if (t == 0) t = 1;

    int ret = pos ^ t;

    ret = ret - bias;
    if (ret < 0) ret += n;
    
    return ret;
}

int alt(int pos, int fp) // get alternate position
{
    static const int len[4] = {32768, 64, 64, 64};

    int fp_hash = fp * 0x5bd1e995;

    //int bias = fp_hash & (max_2_power - 1);

    int bias = fp_hash & (len[0] - 1);
    int segment_length = len[fp & 3];

    //int segment_length = 64;
    //if ((fp & 3) == 0) segment_length = 2048;

    int t = (fp_hash & (segment_length - 1));
    if (t == 0) t = 1;

    pos = pos + bias;
    if (pos >= n) pos -= n;
    
    //int segment_length = 1 << find_highest_bit(pos ^ n);
    //int segment_length = lower_power_of_two((pos ^ n));
    int ret = pos ^ t;

    ret = ret - bias;
    if (ret < 0) ret += n;
    
    return ret;
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

int main()
{
    printf("simulate for cuckoo/vacuum lookup\n");
    const int fp_len = 12;

    // first check !

    /*
    for (int i = 0; i < n; i++)
    {
        int x = i & ((1 << fp_len) - 1);
        if (x == 0) x = 1;
        int t1 = benchmark(i, x);
        int t2 = alt(i, x);
        assert(t1 == t2);
    }
    */

    //srand(time(0));
    srand(233);
    int *T = (int *)malloc(n * 4);
    for (int i = 0; i < n; i += 23) T[i] = rand();
    /*
    int mod = 6789;

    int pw = 1;
    for (; pw <= mod; pw <<= 1);
    pw--;

    deln(pw);
    */

    int sb = 0;
    auto start = chrono::steady_clock::now();
    for (int j = 0; j < n; j++)
    {
        int i = (j * 0x7b9c8a33) & (n - 1); // a random number as the primary bucket position

        int x = i & ((1 << fp_len) - 1); // fingerprint 
        if (x == 0) x = 1; 

        sb += T[i]; // the first memory access

        int t = benchmark(i, x); // alternate function 
        sb += T[t]; // the alternate memory access
    }
    auto end = chrono::steady_clock::now();
    double cost = time_cost(start, end);
    double mops = n / 1000000.0 / cost;

    printf("ignore the value...%d\n", sb);
    printf("benchmark cuckoo: \n");
    printf("time = %.2f, mops = %.2f\n", cost, mops);

    int *T1 = (int *)malloc(n * 4);
    for (int i = 0; i < n; i += 23) T1[i] = rand();
    sb = 0;
    start = chrono::steady_clock::now();
    for (int j = 0; j < n; j++)
    {
        int i = (j * 0x7b9c8a33) % n ; // a random number as the primary bucket position
        if (i < 0) i += n;

        int x = i & ((1 << fp_len) - 1); 
        if (x == 0) x = 1; // fingerprint

        sb += T1[i];  // the first memory access

        int t = alt(i, x); // our alternate function

        if (T1[t] & 1) sb += T1[t]; // the second memory access

        //printf("i = %d, alt = %d\n", i, t);
    }

    end = chrono::steady_clock::now();
    cost = time_cost(start, end);
    mops = n / 1000000.0 / cost;

    printf("ignore the value...%d\n", sb);
    printf("new vacuum: \n");
    printf("time = %.2f, mops = %.2f\n", cost, mops);

    free(T);
    free(T1);

    return 0;
}
