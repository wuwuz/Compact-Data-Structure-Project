#ifndef FILTER_H
#define FILTER_H

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

class Filter
{
    public : 

    void initialization(int _n, int _m, int _max_kick_steps);
    int insert(int ele);
    bool lookup(int ele);
    double get_load_factor();

    protected : 

    int n; // number of buckets
    int m; // number of slots per bucket
    int filled_cell;
    unsigned char *T; // table
    int max_kick_steps;

    int position_hash(int ele); // hash to range [0, n - 1]
    uint8_t fingerprint(int ele); // 32-bit to 8-bit fingerprint
    uint8_t get_item(int pos, int rk); // get cell value (one fingerprint) from tabel[pos][rk]
    void set_item(int pos, int rk, uint8_t fp); // set tabel[pos][rk] as fp

    virtual int alternate(int pos, uint8_t fp){ return 0; }; // get alternate position
    int insert_to_bucket(int pos, uint8_t fp); // insert one fingerprint to bucket [pos] 
    bool lookup_in_bucket(int pos, uint8_t fp); // lookup one fingerprint in  bucket [pos]
} ; 

int Filter::position_hash(int ele)
{
    return (ele % n + n) % n;
}

void Filter::initialization(int _n, int _m, int _max_kick_steps)
{
    n = _n;
    m = _m;
    max_kick_steps = _max_kick_steps;
    filled_cell = 0;

    T = (uint8_t*) calloc(n * m, sizeof(uint8_t));
} 

uint8_t Filter::fingerprint(int ele)
{
    uint8_t h = HashUtil::BobHash32(&ele, 4) % 255 + 1;
    return h;
}

uint8_t Filter::get_item(int pos, int rk)
{
    return (uint8_t) T[pos * m + rk];
}

void Filter::set_item(int pos, int rk, uint8_t fp)
{
    T[pos * m + rk] = (uint8_t) fp;
}

int Filter::insert_to_bucket(int pos, uint8_t fp)
{

    // if success return 0
    // if fail return 1

    for (int i = 0; i < m; i++) 
        if (get_item(pos, i) == 0)
        {
            set_item(pos, i, fp);
            return 0;
        } else 
        {
            //uint8_t t = get_item(pos, i);
        }

    return 1;
}

int Filter::insert(int ele)
{

    // If insert success return 0
    // If insert fail return 1
    //

    uint8_t fp = fingerprint(ele);
    int cur1 = position_hash(ele);
    int cur2 = alternate(cur1, fp);
    if (alternate(cur2, fp) != cur1)
    {
        deln(ele);
        deln((uint8_t)fp);
    }

    if (insert_to_bucket(cur1, fp) == 0) {filled_cell++; return 0;}
    if (insert_to_bucket(cur2, fp) == 0) {filled_cell++; return 0;}

    //randomly choose one bucket's elements to kick
	int cur = (rand() & 1) ? cur1 : cur2;
    int rk = rand() % m;

    //get those item
    uint8_t tmp_fp = get_item(cur, rk);
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

bool Filter::lookup_in_bucket(int pos, uint8_t fp)
{
    // If lookup success return true
    // If lookup fail return false

    for (int i = 0; i < m; i++)
    {
        if (get_item(pos, i) == (uint8_t) fp) 
            return true;
    }
    return false;
}

bool Filter::lookup(int ele)
{

    // If ele is positive, return true
    // negative -- return false

    uint8_t fp = fingerprint(ele);
    int pos1 = position_hash(ele);
    int pos2 = alternate(pos1, fp);

    bool ok1 = lookup_in_bucket(pos1, fp);
    bool ok2 = lookup_in_bucket(pos2, fp);

    return ok1 | ok2;
}

double Filter::get_load_factor()
{
    return filled_cell * 1.0 / n / m;
}

class CuckooFilter : public Filter
{
    private : 
    int alternate(int pos, uint8_t fp) // get alternate position
    {
        return pos ^ position_hash(HashUtil::MurmurHash32(fp));
    }
}; 

class MyFilter : public Filter
{
    private : 

    int alternate(int pos, uint8_t fp)
    {
        // My cuckoo filter -- plus or minus
        int bias = (n / 256) * fp;  // bias -- avoid aggregate /////////////////////
        pos = (pos + bias) % n;
        int delta = HashUtil::MurmurHash32(fp) % (n / 2);
        delta = (delta == 0) ? 1 : delta;

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
					//delta >>= 1;
					//if (delta == 0) delta = 1;
                    delta = HashUtil::MurmurHash32(delta) % nn;
                    delta = (delta == 0) ? 1 : delta;
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
