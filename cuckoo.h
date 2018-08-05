#ifndef CUCKOO_H
#define CUCKOO_H

#include "hashutil.h"
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <algorithm>
using namespace std;

#define memcle(a) memset(a, 0, sizeof(a))
#define sqr(a) ((a) * (a))
#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a) cerr << #a << " = " << a << endl

template <typename fp_t, int fp_len>
class Filter
{
    public : 

    int n; // number of buckets
    int m; // number of slots per bucket
    void init(int _n, int _m, int _max_kick_steps);
	void clear();
    int insert(int ele);
    bool lookup(int ele);
    double get_load_factor();

	~Filter() { free(T); }

    protected : 

    int filled_cell;
	fp_t *T; // hash table
    int max_kick_steps;

    int position_hash(int ele); // hash to range [0, n - 1]
    fp_t fingerprint(int ele); // 32-bit to 'fp_len'-bit fingerprint
    fp_t get_item(int pos, int rk); // get cell value (one fingerprint) from tabel[pos][rk]
    void set_item(int pos, int rk, fp_t fp); // set tabel[pos][rk] as fp

    virtual int alternate(int pos, fp_t fp){ return 0; }; // get alternate position
    int insert_to_bucket(int pos, fp_t fp); // insert one fingerprint to bucket [pos] 
    bool lookup_in_bucket(int pos, fp_t fp); // lookup one fingerprint in  bucket [pos]
} ; 

template <typename fp_t, int fp_len>
int Filter<fp_t, fp_len>::position_hash(int ele)
{
    return (ele % n + n) % n;
}

template <typename fp_t, int fp_len>
void Filter<fp_t, fp_len>::init(int _n, int _m, int _step)
{
    n = _n;
    m = _m;
    max_kick_steps = _step;
    filled_cell = 0;

    T = (fp_t*) calloc(n * m, sizeof(fp_t));
}

template <typename fp_t, int fp_len>
void Filter<fp_t, fp_len>::clear()
{
	filled_cell = 0;
	memset(T, 0, sizeof(fp_t) * (n * m));
}

template <typename fp_t, int fp_len>
fp_t Filter<fp_t, fp_len>::fingerprint(int ele)
{
    fp_t h = HashUtil::BobHash32(&ele, 4) % ((1ull << fp_len) -1) + 1;
    return h;
}

template <typename fp_t, int fp_len>
fp_t Filter<fp_t, fp_len>::get_item(int pos, int rk)
{
    return (fp_t) T[pos * m + rk];
}

template <typename fp_t, int fp_len>
void Filter<fp_t, fp_len>::set_item(int pos, int rk, fp_t fp)
{
    T[pos * m + rk] = (fp_t) fp;
}

template <typename fp_t, int fp_len>
int Filter<fp_t, fp_len>::insert_to_bucket(int pos, fp_t fp)
{

    // if success return 0
    // if fail return 1

    for (int i = 0; i < m; i++) 
        if (get_item(pos, i) == 0)
        {
            set_item(pos, i, fp);
            return 0;
        }

    return 1;
}

template <typename fp_t, int fp_len>
int Filter<fp_t, fp_len>::insert(int ele)
{

    // If insert success return 0
    // If insert fail return 1

    fp_t fp = fingerprint(ele);
    int cur1 = position_hash(ele);
    int cur2 = alternate(cur1, fp);
    if (alternate(cur2, fp) != cur1)
    {
        deln(ele);
        deln((fp_t)fp);
    }

    if (insert_to_bucket(cur1, fp) == 0) {filled_cell++; return 0;}
    if (insert_to_bucket(cur2, fp) == 0) {filled_cell++; return 0;}

    //randomly choose one bucket's elements to kick
	int cur = (rand() & 1) ? cur1 : cur2;
    int rk = rand() % m;

    //get those item
    fp_t tmp_fp = get_item(cur, rk);
    set_item(cur, rk, fp);

    int alt = alternate(cur, tmp_fp);
    
    for (int i = 0; i < max_kick_steps; i++)
    {
        if (insert_to_bucket(alt, tmp_fp) == 0) {filled_cell++; return 0;}
        rk = rand() % m;
        fp = get_item(alt, rk);
        set_item(alt, rk, tmp_fp);

        tmp_fp = fp;
        alt = alternate(alt, tmp_fp);
    }

    return 1;
}

template <typename fp_t, int fp_len>
bool Filter<fp_t, fp_len>::lookup_in_bucket(int pos, fp_t fp)
{
    // If lookup success return true
    // If lookup fail return false

    for (int i = 0; i < m; i++)
    {
        if (get_item(pos, i) == (fp_t) fp) 
            return true;
    }
    return false;
}

template <typename fp_t, int fp_len>
bool Filter<fp_t, fp_len>::lookup(int ele)
{

    // If ele is positive, return true
    // negative -- return false

	fp_t fp = fingerprint(ele);
    int pos1 = position_hash(ele);
    int pos2 = alternate(pos1, fp);

    bool ok1 = lookup_in_bucket(pos1, fp);
    bool ok2 = lookup_in_bucket(pos2, fp);

    return ok1 | ok2;
}

template <typename fp_t, int fp_len>
double Filter<fp_t, fp_len>::get_load_factor()
{
    return filled_cell * 1.0 / n / m;
}

template <typename fp_t, int fp_len>
class CuckooFilter : public Filter<fp_t, fp_len>
{
    private : 
    int alternate(int pos, fp_t fp) // get alternate position
    {
        return pos ^ this->position_hash(HashUtil::MurmurHash32(fp));
    }
}; 

template <typename fp_t, int fp_len>
class MyFilter : public Filter<fp_t, fp_len>
{
    private : 

    int alternate(int pos, fp_t fp)
    {
        // My cuckoo filter -- plus or minus
		int n = this->n;
        int bias = (n / (1 << fp_len)) * fp;  // bias -- avoid aggregate
		pos = (pos + bias) % n;
        int delta = HashUtil::MurmurHash32(fp) % (n-1) + 1;

        int sum_bound = 0;
        int ret = 0;

        int nn = n;

        for (; ;)
        {
            if (((pos / delta) & 1) == 0) 
            {
                if (pos + delta < nn)
                {
                    ret = pos + delta;
                    break;
                } else
                {
                    int lb = max((pos / delta) * delta, nn - delta);
                    int ub = min((pos / delta) * delta + delta - 1, nn - 1); // [lb, ub] is the interval which can't determine pair
                    sum_bound += lb;
                    nn = ub - lb + 1;
                    delta = HashUtil::MurmurHash32(delta) % (nn-1) + 1;
                    pos = pos - lb;
                }
            }
            else
            {
                ret = pos - delta; // case : minus delta -- must be legal !
                break;
            }
        }

        ret = ret + sum_bound;
        ret = (ret - bias + n) % n;
        return ret;
    }
};

#endif
