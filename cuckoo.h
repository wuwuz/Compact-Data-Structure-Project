#ifndef CUCKOO_H
#define CUCKOO_H

#include "hashutil.h"
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <random>
using namespace std;

#define memcle(a) memset(a, 0, sizeof(a))
#define sqr(a) ((a) * (a))
#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a) cerr << #a << " = " << a << endl
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ROUNDDOWN(a, b) ((a) - ((a) % (b)))
#define ROUNDUP(a, b) ROUNDDOWN((a) + (b - 1), b)

inline int find_the_highest_bit(int v)
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



template <typename fp_t, int fp_len>
class Filter
{
    public : 

    int n; // number of buckets
    int m; // number of slots per bucket
    int memory_consumption;
    virtual void init(int _n, int _m, int _max_kick_steps) = 0;
	virtual void clear() = 0;
    virtual int insert(int ele) = 0;
    virtual bool lookup(int ele) = 0;
    int position_hash(int ele); // hash to range [0, n - 1]
    virtual double get_load_factor() = 0;
    virtual double get_full_bucket_factor() = 0;
    virtual void debug_test() = 0;
};

template <typename fp_t, int fp_len>
int Filter<fp_t, fp_len>::position_hash(int ele)
{
    return (ele % n + n) % n;
}

template <typename fp_t, int fp_len>
class BloomFilter : public Filter<fp_t, fp_len>
{
    public : 

    int shift;
    fp_t *T;
	~BloomFilter() { free(T); }

    void init(int _n, int _m, int _max_kick_steps = 0)
    {
        shift = 0;
        int t = fp_len;
        for (; t > 1; t >>= 1) shift++;

        this -> n = _n * _m * fp_len;
        this -> memory_consumption = _n * _m * sizeof(fp_t);
        T = (fp_t*) calloc(_n * _m, sizeof(fp_t)); // how many bytes
    }
    void clear()
    {
	    memset(T, 0, sizeof(fp_t) * (this -> n / fp_len));
    }
    void set_item(int pos)
    {
        T[pos >> shift] |= (1 << (pos & ((1 << shift) - 1)));
    }
    bool get_item(int pos)
    {
        return (T[pos >> shift] & (1 << (pos & ((1 << shift) - 1)))) > 0;
    }
    int insert(int ele)
    {
        int h1 = HashUtil::MurmurHash32(ele);
        int h2 = HashUtil::MurmurHash32(h1);
        int h3 = HashUtil::MurmurHash32(h2);

        set_item(this->position_hash(h1));
        set_item(this->position_hash(h2));
        set_item(this->position_hash(h3));
        return 0;
    }
    bool lookup(int ele)
    {
        int h1 = HashUtil::MurmurHash32(ele);
        int h2 = HashUtil::MurmurHash32(h1);
        int h3 = HashUtil::MurmurHash32(h2);
        
        return get_item(this->position_hash(h1)) && get_item(this->position_hash(h2)) && get_item(this->position_hash(h3));
    }
    double get_load_factor(){return 0;}
    double get_full_bucket_factor(){return 0;}
    void debug_test() {}
};

template <typename fp_t, int fp_len>
class CuckooFilter : public Filter<fp_t, fp_len>
{
    public : 

    int max_2_power;
    virtual void init(int _n, int _m, int _max_kick_steps);
	void clear();
    int insert(int ele);
    bool lookup(int ele);
    double get_load_factor();

    fp_t *T;
	~CuckooFilter() { free(T); }
    
    int filled_cell;
    int max_kick_steps;

    fp_t fingerprint(int ele); // 32-bit to 'fp_len'-bit fingerprint
    virtual fp_t get_item(int pos, int rk); // get cell value (one fingerprint) from tabel[pos][rk]
    virtual void set_item(int pos, int rk, fp_t fp); // set tabel[pos][rk] as fp

    virtual int alternate(int pos, fp_t fp) = 0; // get alternate position
    virtual int insert_to_bucket(int pos, fp_t fp); // insert one fingerprint to bucket [pos] 
    virtual int lookup_in_bucket(int pos, fp_t fp); // lookup one fingerprint in  bucket [pos]
    virtual void debug_test() {}
} ; 


template <typename fp_t, int fp_len>
void CuckooFilter<fp_t, fp_len>::init(int _n, int _m, int _step)
{
    this -> n = _n;
    this -> m = _m;
    this -> max_kick_steps = _step;
    this -> filled_cell = 0;
    this -> memory_consumption = int(_n * _m * (fp_len) * 1.0 / 8);

    max_2_power = 1;
    for (; max_2_power * 2 < _n; ) max_2_power <<= 1;

    this -> T = (fp_t*) calloc(this -> memory_consumption, sizeof(char));
}

template <typename fp_t, int fp_len>
void CuckooFilter<fp_t, fp_len>::clear()
{
	this -> filled_cell = 0;
	memset(this -> T, 0, this -> memory_consumption);
}

template <typename fp_t, int fp_len>
fp_t CuckooFilter<fp_t, fp_len>::fingerprint(int ele)
{
    fp_t h = HashUtil::BobHash32(&ele, 4) % ((1ull << fp_len) -1) + 1;
    return h;
}

template <typename fp_t, int fp_len>
fp_t CuckooFilter<fp_t, fp_len>::get_item(int pos, int rk)
{
    return (fp_t) this -> T[pos * this -> m + rk];
}

template <typename fp_t, int fp_len>
void CuckooFilter<fp_t, fp_len>::set_item(int pos, int rk, fp_t fp)
{
    this -> T[pos * this -> m + rk] = (fp_t) fp;
}

template <typename fp_t, int fp_len>
int CuckooFilter<fp_t, fp_len>::insert_to_bucket(int pos, fp_t fp)
{

    // if success return 0
    // if fail return 1

    for (int i = 0; i < this -> m; i++) 
        if (get_item(pos, i) == 0)
        {
            set_item(pos, i, fp);
            return 0;
        }

    return 1;
}

template <typename fp_t, int fp_len>
int CuckooFilter<fp_t, fp_len>::insert(int ele)
{

    // If insert success return 0
    // If insert fail return 1

    fp_t fp = fingerprint(ele);
    int cur1 = this -> position_hash(ele);
    int cur2 = alternate(cur1, fp);
    if (alternate(cur2, fp) != cur1)
    {
        deln(ele);
        deln(fp);
        printf("cur1 = %d alt cur1 = %d\n", cur1, alternate(cur1, fp));
        printf("cur2 = %d alt cur2 = %d\n", cur2, alternate(cur2, fp));
    }

    if (insert_to_bucket(cur1, fp) == 0) {filled_cell++; return 0;}
    if (insert_to_bucket(cur2, fp) == 0) {filled_cell++; return 0;}

    //randomly choose one bucket's elements to kick
	int cur = (rand() & 1) ? cur1 : cur2;
    int rk = rand() % this -> m;

    //get those item
    fp_t tmp_fp = get_item(cur, rk);
    set_item(cur, rk, fp);

    int alt = alternate(cur, tmp_fp);
    
    for (int i = 0; i < this -> max_kick_steps; i++)
    {
        if (insert_to_bucket(alt, tmp_fp) == 0) {filled_cell++; return 0;}
        rk = rand() % this -> m;
        fp = get_item(alt, rk);
        set_item(alt, rk, tmp_fp);

        tmp_fp = fp;
        alt = alternate(alt, tmp_fp);
    }

    return 1;
}

template <typename fp_t, int fp_len>
int CuckooFilter<fp_t, fp_len>::lookup_in_bucket(int pos, fp_t fp)
{
    // If lookup success return 1
    // If lookup fail and the bucket is full return 2
    // If lookup fail and the bucket is not full return 3

    int isFull = 1;
    for (int i = 0; i < this -> m; i++)
    {
        int t = get_item(pos, i);
        if (t == (fp_t) fp) 
            return 1;
        isFull &= (t != 0);
    }
    return (isFull) ? 2 : 3;
}

template <typename fp_t, int fp_len>
bool CuckooFilter<fp_t, fp_len>::lookup(int ele)
{

    // If ele is positive, return true
    // negative -- return false

	fp_t fp = fingerprint(ele);
    int pos1 = this -> position_hash(ele);
    int pos2 = alternate(pos1, fp);

    int ok1 = lookup_in_bucket(pos1, fp);

    if (ok1 == 1) return true;
    if (ok1 == 3) return false;

    int ok2 = lookup_in_bucket(pos2, fp);

    return ok2 == 1;
}

template <typename fp_t, int fp_len>
double CuckooFilter<fp_t, fp_len>::get_load_factor()
{
    return filled_cell * 1.0 / this -> n / this -> m;
}

template <typename fp_t, int fp_len>
class SemiSortCuckooFilter : public Filter<fp_t, fp_len>
{
    public : 

    int max_2_power;
    virtual void init(int _n, int _m, int _max_kick_steps);
	void clear();
    int insert(int ele);
    bool lookup(int ele);
    double get_load_factor();
    double get_full_bucket_factor();

    bool debug_flag = false;
    bool balance = true;
    uint32_t *T;
    uint32_t encode_table[1 << 16];
    uint32_t decode_table[1 << 16];

	~SemiSortCuckooFilter() { free(T); }
    
    int filled_cell;
    int full_bucket;
    int max_kick_steps;

    fp_t fingerprint(int ele); // 32-bit to 'fp_len'-bit fingerprint

    //interface for semi-sorted bucket
    void get_bucket(int pos, fp_t *store);
    void set_bucket(int pos, fp_t *sotre);
    void test_bucket();
    void make_balance();
    int high_bit(fp_t fp);
    int low_bit(fp_t fp);

    virtual int alternate(int pos, fp_t fp) = 0; // get alternate position
    virtual int insert_to_bucket(fp_t *store, fp_t fp); // insert one fingerprint to bucket [pos] 
    virtual int lookup_in_bucket(int pos, fp_t fp); // lookup one fingerprint in  bucket [pos]
    void debug_test()
    {
        //debug_flag = true;
        //static int in_deg[1000][1000];
        //memcle(in_deg);

        /*
        mt19937 rd(123);
        for (int i = 1; i <= 400; i++)
        {
            int pos = this -> position_hash(rd());
            fp_t store[8];
            get_bucket(pos, store);

            int t = 0;
            for (int j = 0; j < this -> m; j++)
                t += store[j] != 0;
            if (t != this -> m)
                printf("pos = %d, cnt = %d\n", pos, t);
        }
        */

        /*
        static int in_deg[1000][1000];
        memset(in_deg, 1, sizeof(in_deg));
        for (int i = 0; i < this -> n; i++) in_deg[i][i] = 0;

        //printf("non-full index : \n");
        for (int i = 0; i < this -> n; i++)
        {
            fp_t store[8];
            get_bucket(i, store);
            int t = 0;
            for (int j = 0; j < this -> m; j++)
                if (store[j] != 0)
                {
                    int alt = this -> alternate(i, store[j]);
                    in_deg[i][alt] = 1;
                    t += store[j] != 0;
                }
            if (t != this -> m)
                printf("pos = %d, cnt = %d\n", i, t);
        }
        puts("");

        for (int k = 0; k < this -> n; k++)
            for (int i = 0; i < this -> n; i++)
                if (i != k)
                    for (int j = 0; j < this -> n; j++)
                        if (j != k && i != j)
                            in_deg[i][j] = MIN(in_deg[i][j], in_deg[i][k] + in_deg[k][j]);

            */

        /*

        for (int i = 0; i < this -> n; i++)
            for (int j = 0; j < (1 << fp_len); j++)
            {
                int t = this -> alternate(i, j);
                //printf("t = %d, i = %d\n", t, i);
                in_deg[i][t]++;
            }
        
        for (int i = 0; i < this -> n; i++)
        {
            for (int j = 0; j < this -> n; j++)
                printf("%d ", in_deg[i][j]);
            puts("");
        }
        */
    }
} ; 


template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::init(int _n, int _m, int _step)
{
    this -> n = _n;
    this -> m = _m;
    this -> max_kick_steps = _step;
    this -> filled_cell = 0;
    this -> full_bucket = 0;

    int how_many_bit = this -> n * this -> m * (fp_len - 1);

    this -> memory_consumption = ROUNDUP(ROUNDUP(how_many_bit, 8) / 8 + 8, 8); // how many bytes !

    max_2_power = 1;
    for (; max_2_power * 2 < _n; ) max_2_power <<= 1;

    //this -> T = (fp_t*) calloc(this -> n * this -> m, sizeof(fp_t));
    this -> T = (uint32_t *) calloc(this -> memory_consumption, sizeof(char));

    if (this -> m == 4)
    {
        int index = 0;
        /*
        for (int i = 0; i < 16; i++)
            for (int j = 0; j < ((i == 0) ? 1 : i); j++)
                for (int k = 0; k < ((j == 0) ? 1 : j); k++)
                    for (int l = 0; l < ((k == 0) ? 1 : k); l++)
                    {
                        int plain_bit = (i << 12) + (j << 8) + (k << 4) + l;
                        encode_table[plain_bit] = index;
                        decode_table[index] = plain_bit;
                        ++index;
                    }
        */
        for (int i = 0; i < 16; i++)
            for (int j = 0; j < ((i == 0) ? 1 : i + 1); j++)
                for (int k = 0; k < ((j == 0) ? 1 : j + 1); k++)
                    for (int l = 0; l < ((k == 0) ? 1 : k + 1); l++)
                    {
                        int plain_bit = (i << 12) + (j << 8) + (k << 4) + l;
                        encode_table[plain_bit] = index;
                        decode_table[index] = plain_bit;
                        ++index;
                    }
        deln(index);
    }
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::clear()
{
	this -> filled_cell = 0;
	//memset(this -> T, 0, sizeof(fp_t) * (this -> n * this -> m));
    memset(this -> T, 0, this -> memory_consumption);
}

template <typename fp_t, int fp_len>
fp_t SemiSortCuckooFilter<fp_t, fp_len>::fingerprint(int ele)
{
    fp_t h = HashUtil::BobHash32(&ele, 4) % ((1ull << fp_len) -1) + 1;
    return h;
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::get_bucket(int pos, fp_t *store)
{
    // Default : 
    //
    // Little Endian Store
    // Store by uint32_t


    // 1. read the endcoded bits from memory

    int bucket_length = (fp_len - 1) * this -> m;
    int start_bit_pos = pos * bucket_length;
    int end_bit_pos = start_bit_pos + bucket_length - 1;
    uint64_t result = 0;
    int readed_bit_count = 0;

    for (int i = start_bit_pos; i < end_bit_pos; ) // i : bit position
    {
        // in this 32-bit unit : 
        uint32_t unit = ((uint32_t *)T)[ROUNDDOWN(i, 32) / 32];
        int reading_lower_bound = MAX(start_bit_pos, ROUNDDOWN(i, 32)) & 31;
        int reading_upper_bound = MIN(end_bit_pos, ROUNDDOWN(i, 32) + 31) & 31;

        uint64_t reading_result = ((uint64_t)unit & ((1LL << (reading_upper_bound + 1)) - 1)) >> reading_lower_bound;

        result = (reading_result << readed_bit_count) + result;
        readed_bit_count += (reading_upper_bound - reading_lower_bound + 1);
        i += (reading_upper_bound - reading_lower_bound + 1);
    }


    // 2. read the 4 elements from the encoded bits
    // We use 12 bits to store the 16 most significant bits for the items in bucket, 4 bits per item
    // the low bits are stored in the remaining bits
    //
    // For example, 8 bits per item , require 28 bits to store: 
    //
    // Original : 
    //
    // hhhh llll
    // hhhh llll
    // hhhh llll
    // hhhh llll
    //
    // encoded : 
    //
    //
    // 0 - 11                       12 - 15    16 - 19  20-23   24 - 27
    // HHHHHHHHHHHH                 llll       llll     llll    llll
    //  encoded high bit(12 bits)   item 0     item 1   item 2  item 3
    //

    for (int i = this -> m - 1; i >= 0; i--)
    {
        store[i] = result & ((1 << (fp_len - 4)) - 1);
        result >>= (fp_len - 4);
    }

    int decode_result = decode_table[result];
    //printf("high_code = %x\n", decode_result);

    for (int i = this -> m - 1; i >= 0; i--)
    {
        store[i] = ((decode_result & ((1 << 4) - 1)) << (fp_len - 4)) + store[i];
        decode_result >>= 4;
    }
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::set_bucket(int pos, fp_t *store)
{
    // 0. sort store ! descendant order >>>>>>

    for (int i = 0; i < this -> m; i++)
        for (int j = i + 1; j < this -> m; j++)
            if (store[j] > store[i])
                swap(store[i], store[j]);
    
    /*
    for (int i = 0; i < this -> m - 1; i++)
        if (store[i] != 0 && store[i + 1] == store[i])
            printf("same ? ");
    */
    // 1. compute the encode 

    uint64_t high_bit = 0;
    uint64_t low_bit = 0;

    for (int i = 0; i < this -> m; i++)
    {
        low_bit = (low_bit << (fp_len - 4)) + (store[i] & ((1 << (fp_len - 4)) - 1));
        high_bit = (high_bit << 4) + ((store[i] >> (fp_len - 4)) & ((1 << 4) - 1));
    }

    // 2. store into memory
    uint64_t high_encode = encode_table[high_bit];
    //printf("high_bit = %x\n", high_bit);
    //printf("high_encode = %x\n", high_encode);
    uint64_t all_encode = (high_encode << (this -> m * (fp_len - 4))) + low_bit;

    int bucket_length = (fp_len - 1) * this -> m;
    int start_bit_pos = pos * bucket_length;
    int end_bit_pos = start_bit_pos + bucket_length - 1;

    for (int i = start_bit_pos; i < end_bit_pos; ) // i : bit position
    {
        // in this 32-bit unit : 
        uint32_t unit = ((uint32_t *)T)[ROUNDDOWN(i, 32) / 32];
        int writing_lower_bound = MAX(start_bit_pos, ROUNDDOWN(i, 32)) & 31;
        int writing_upper_bound = MIN(end_bit_pos, ROUNDDOWN(i, 32) + 31) & 31;

        uint64_t tmp = (unit & ((1LL << writing_lower_bound) - 1));
        tmp += ((all_encode & ((1LL << (writing_upper_bound - writing_lower_bound + 1)) - 1)) << writing_lower_bound);
        //tmp += (unit & (~((1 << (writing_upper_bound + 1)) - 1)));
        uint32_t sb =  ((uint64_t)unit >> (writing_upper_bound + 1)) << (writing_upper_bound + 1);
        tmp += sb;

        ((uint32_t *)T)[ROUNDDOWN(i, 32) / 32] = tmp;

        i += (writing_upper_bound - writing_lower_bound + 1);
        all_encode >>= (writing_upper_bound - writing_lower_bound + 1);
    }

}

/*
template <typename fp_t, int fp_len>
fp_t CuckooFilter<fp_t, fp_len>::get_item(int pos, int rk)
{
    return (fp_t) this -> T[pos * this -> m + rk];
}

template <typename fp_t, int fp_len>
void CuckooFilter<fp_t, fp_len>::set_item(int pos, int rk, fp_t fp)
{
    this -> T[pos * this -> m + rk] = (fp_t) fp;
}
*/

template <typename fp_t, int fp_len>
inline int SemiSortCuckooFilter<fp_t, fp_len>::high_bit(fp_t fp)
{
    return (fp >> (fp_len - 4)) & ((1 << 4) - 1);
}

template <typename fp_t, int fp_len>
inline int SemiSortCuckooFilter<fp_t, fp_len>::low_bit(fp_t fp)
{
    return fp & ((1 << (fp_len - 4)) - 1);
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::insert_to_bucket(fp_t *store, fp_t fp)
{

    // if success return 0
    // if find the same : return 1 + position
    // if full : return 1 + 4

    //get_bucket(pos, store);

    for (int i = 0; i < this -> m; i++) 
    {
        if (store[i] == 0)
        {
            if (i == this -> m - 1) // full
                full_bucket++;
            store[i] = fp;
            //set_bucket(pos, store);
            return 0;
        } 
        /*
        else
        // optimization for one bit
        if (high_bit(store[i]) == high_bit(fp))
            return 1 + i;
        */
    }

    return 1 + 4;
}

/*
template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::insert(int ele)
{

    // If insert success return 0
    // If insert fail return 1

    fp_t fp = fingerprint(ele);
    int cur1 = this -> position_hash(ele);
    int cur2 = alternate(cur1, fp);
    if (alternate(cur2, fp) != cur1)
    {
        deln(ele);
        deln(fp);
        printf("cur1 = %d alt cur1 = %d\n", cur1, alternate(cur1, fp));
        printf("cur2 = %d alt cur2 = %d\n", cur2, alternate(cur2, fp));
    }

    fp_t store1[8];
    fp_t store2[8];

    get_bucket(cur1, store1);
    get_bucket(cur2, store2);

    if (insert_to_bucket(store1, fp) == 0) {filled_cell++; set_bucket(cur1, store1); return 0;}
    if (insert_to_bucket(store2, fp) == 0) {filled_cell++; set_bucket(cur2, store2); return 0;}

    //randomly choose one bucket's elements to kick
    int rk = rand() % this -> m;

    //get those item
    int cur;
    fp_t *cur_store;

    if (rand() & 1)
        cur = cur1, cur_store = store1;
    else 
        cur = cur2, cur_store = store2;

    fp_t tmp_fp = cur_store[rk];
    cur_store[rk] = fp;
    set_bucket(cur, cur_store);

    int alt = alternate(cur, tmp_fp);
    
    for (int i = 0; i < this -> max_kick_steps; i++)
    {
        memset(store1, 0, sizeof(store1));
        get_bucket(alt, store1);
        if (insert_to_bucket(store1, tmp_fp) == 0) 
        {
            filled_cell++; 
            set_bucket(alt, store1);
            return 0;
        }

        rk = rand() % this -> m;
        fp = store1[rk];
        store1[rk] = tmp_fp;
        set_bucket(alt, store1);

        tmp_fp = fp;
        alt = alternate(alt, tmp_fp);
    }

    return 1;
}
*/

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::make_balance()
{
    int success = 0;
    for (int i = 0; i < this -> n; i++)
    {
        fp_t store[8], tmp[8];
        get_bucket(i, store);
        if (store[this -> m - 1] != 0) // full
        {
            for (int j = 0; j < this -> m; j++)
            {
                int alt = this -> alternate(i, store[j]);
                get_bucket(alt, tmp);
                if (tmp[this -> m - 2] == 0) // not full, <= 2
                {
                    tmp[this -> m - 2] = store[j];
                    set_bucket(alt, tmp);
                    store[j] = 0;
                    set_bucket(i, store);
                    full_bucket--;
                    success++;
                    break;
                }
            }
        }
    }
    deln(success);
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::insert(int ele)
{

    if (SemiSortCuckooFilter<fp_t, fp_len>::lookup(ele) == true)
    {
        ++filled_cell;
        return 0;
    }
    // If insert success return 0
    // If insert fail return 1

    /*
    if (balance == false && get_load_factor() > 0.94)
    {
        make_balance();
        balance = true;
    }
    */

    fp_t fp = fingerprint(ele);
    int cur1 = this -> position_hash(ele);
    int cur2 = alternate(cur1, fp);
    if (alternate(cur2, fp) != cur1)
    {
        deln(ele);
        deln(fp);
        printf("cur1 = %d alt cur1 = %d\n", cur1, alternate(cur1, fp));
        printf("cur2 = %d alt cur2 = %d\n", cur2, alternate(cur2, fp));
    }

    fp_t store1[8];
    fp_t store2[8];

    int cnt1 = 0, cnt2 = 0;
    int res1 = 0, res2 = 0;
    get_bucket(cur1, store1);
    for (int i = 0; i < this -> m; i++) cnt1 += (store1[i] != 0);
    get_bucket(cur2, store2);
    for (int i = 0; i < this -> m; i++) cnt2 += (store2[i] != 0);

    if (cnt1 <= cnt2)
    {
        if ((res1 = insert_to_bucket(store1, fp)) == 0) {filled_cell++; set_bucket(cur1, store1); return 0;}
        //if ((res2 = insert_to_bucket(store2, fp)) == 0) {filled_cell++; set_bucket(cur2, store2); return 0;}
    } else
    {
        if ((res2 = insert_to_bucket(store2, fp)) == 0) {filled_cell++; set_bucket(cur2, store2); return 0;}
        //if ((res1 = insert_to_bucket(store1, fp)) == 0) {filled_cell++; set_bucket(cur1, store1); return 0;}
    }

    /*
    //choose one slot to insert fp
    int cur;
    fp_t *cur_store;

    if (rand() & 1)
        cur = cur1, cur_store = store1;
    else 
        cur = cur2, cur_store = store2;
        */

    // use bfs to search one non-full bucket

    static fp_t bucket[2000][4];
    static int bucket_pos[2000];
    int bucket_cnt = 0;

    static int index[3000];
    static int id[3000];
    static int prev[3000];
    static fp_t que_fp[3000];

    // index[i][0] : index in bucket array
    // id[i][1] : index in one bucket
    // prev[i][2] : previous status in bfs array
    // que_fp[i] : the fp queue!

    bool opt = false;

    if (get_load_factor() > 0.95)
        opt = true;

    int l = 0, r = 0, final_cnt = 100, p = -1;

    //BloomFilter<uint8_t, 8> filter;
    //filter.init(this -> max_kick_steps * 4, 4); //  200 * 4 * 8

    /*
    if (res1 < 5)
    {
        //opt for one bit
        for (int i = 0; i < this -> m; i++) bucket[0][i] = store1[i];
        ++r;
        index[r] = 0;
        id[r] = res1 - 1;
        prev[r] = -1;
        que_fp[r] = store1[id[r]];
        filter.insert((cur1 << 2) + id[r]);
    }
    else 
    {
    */
    for (int i = 0; i < this -> m; i++) 
    {
        bucket[0][i] = store1[i];
        ++r;
        index[r] = 0;
        id[r] = i;
        prev[r] = -1;
        que_fp[r] = store1[i];
        //filter.insert((cur1 << 2) + id[r]);
    }

    bucket_pos[0] = cur1;

    /*
    if (res2 < 5)
    {
        // opt for one bit
        for (int i = 0; i < this -> m; i++) bucket[1][i] = store2[i];
        ++r;
        index[r] = 1;
        id[r] = res2 - 1;
        prev[r] = -1;
        que_fp[r] = store2[id[r]];
        filter.insert((cur2 << 2) + id[r]);
    }
    else
    */
    for (int i = 0; i < this -> m; i++) 
    {
        ++r;
        bucket[1][i] = store2[i];
        index[r] = 1;
        id[r] = i;
        prev[r] = -1;
        que_fp[r] = store2[i];
        //filter.insert((cur2 << 2) + id[r]);
    }

    bucket_pos[1] = cur2;

    bucket_cnt = 1;
    bool quit = false;

    for (; l < r && r < max_kick_steps * 4 && quit == false; )
    {
        ++l;
        fp_t cur_fp = que_fp[l];
        int pos = bucket_pos[index[l]];
        //deln(pos);
        int alt = alternate(pos, cur_fp);
        if (cur_fp == 0) continue;
        //if (filter.lookup(alt) == true) continue;

        int bc = ++bucket_cnt;
        get_bucket(alt, bucket[bc]);
        bucket_pos[bc] = alt;
        //filter.insert(alt);
        //

        /*
        bool collide = false;
        for (int i = 0; i < this -> m; i++)
            if (bucket[bc][i] != 0 && high_bit(bucket[bc][i]) == high_bit(cur_fp))
            {
                collide = true;

                if (filter.lookup((alt << 2) + i) == true) break;
                filter.insert((alt << 2) + i);

                ++r;
                index[r] = bc;
                id[r] = i;
                prev[r] = l;
                que_fp[r] = bucket[bc][i];

                break;
            }

        if (collide == true) continue;
        */

        for (int i = 0; i < this -> m; i++) 
        {
            //if (filter.lookup((alt << 2) + i) == true) continue;
            //filter.insert((alt << 2) + i);

            ++r;
            index[r] = bc;
            id[r] = i;
            prev[r] = l;
            que_fp[r] = bucket[bc][i];

            if (bucket[bc][i] == 0)
            {
                // find a empty slot !
                if (final_cnt > i)
                {
                    final_cnt = i;
                    p = r;
                }
                // if opt : find the best
                // if not opt : find the first
                if (opt == false || final_cnt == 0)  
                {
                    quit = true;
                }

                break;
            }
        }
    }

    /*
    if (bucket_cnt >= max_kick_steps)
    {
        deln(bucket_cnt);
        deln(r);
    }
    */
    /*
    assert(bucket_cnt < max_kick_steps);
    assert(r < max_kick_steps * 4);
    */

    if (final_cnt == 100) 
    {
        /*
        printf("fp = %x\n", fp);
        printf("que pos : \n");
        for (int i = 0; i < bucket_cnt; i++)
            printf("%d ", bucket_pos[i]);
        puts("");
        */
        return 1; // insert fail
    }

    for (; p > 0; p = prev[p])
    {
        if (prev[p] == -1)
            bucket[index[p]][id[p]] = fp; // move the original fp to the bucket
        else
            bucket[index[p]][id[p]] = que_fp[prev[p]]; // move the previous fp into current bucket !

        set_bucket(bucket_pos[index[p]], bucket[index[p]]);
    }

    filled_cell++;
    if (final_cnt == this -> m - 1)
        full_bucket++;

    return 0; // succeed !
    /*
    puts("");
    */
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::lookup_in_bucket(int pos, fp_t fp)
{
    // If lookup success return 1
    // If lookup fail and the bucket is full return 2
    // If lookup fail and the bucket is not full return 3

    fp_t store[8];
    get_bucket(pos, store);

    int isFull = 1;
    for (int i = 0; i < this -> m; i++)
    {
        fp_t t = store[i];
        if (t == fp) 
            return 1;
        isFull &= (t != 0);
    }
    return (isFull) ? 2 : 3;
}

template <typename fp_t, int fp_len>
bool SemiSortCuckooFilter<fp_t, fp_len>::lookup(int ele)
{

    // If ele is positive, return true
    // negative -- return false

	fp_t fp = fingerprint(ele);
    int pos1 = this -> position_hash(ele);
    int pos2 = alternate(pos1, fp);

    int ok1 = lookup_in_bucket(pos1, fp);

    if (ok1 == 1) return true;
    //if (ok1 == 3) return false;

    int ok2 = lookup_in_bucket(pos2, fp);

    return ok2 == 1;
}

template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_load_factor()
{
    return filled_cell * 1.0 / this -> n / this -> m;
}

template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_full_bucket_factor()
{
    return full_bucket * 1.0 / this -> n;
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::test_bucket()
{
    for (int i = 0; i < this -> n; i++)
    {
        fp_t store[8];
        store[0] = rand() % (1 << fp_len);
        store[1] = rand() % (1 << fp_len);
        store[2] = rand() % (1 << fp_len);
        store[3] = rand() % (1 << fp_len);

        for (int j = 0; j < this -> m; j++)
            for (int k = j + 1; k < this -> m; k++)
                if (store[j] == store[k])
                    store[k] = 0;

        set_bucket(i, store);

        fp_t tmp_store[8];
        get_bucket(i, tmp_store);
        for (int j = 0; j < 4; j++)
        {
            if (tmp_store[j] != store[j])
                printf("i = %d, j = %d\n", i, j);
            assert(tmp_store[j] == store[j]);
        }
    }
}

template <typename fp_t, int fp_len>
class MortonFilter : public Filter<fp_t, fp_len>
{
    public : 

    int max_2_power;
    virtual void init(int _n, int _m, int _max_kick_steps);
	void clear();
    int insert(int ele);
    bool lookup(int ele);
    double get_load_factor();

    uint8_t *T;

	~MortonFilter() { free(T); }
    
    int filled_cell;
    int max_kick_steps;
    int block_number;

    fp_t fingerprint(int ele); // 32-bit to 'fp_len'-bit fingerprint

    //interface for semi-sorted bucket
    void get_bucket(int pos, fp_t *store);
    void set_bucket(int pos, fp_t *sotre);
    void test_bucket();
    void print_block(int pos);
    int get_count(int block, int index);
    void set_count(int block, int index, int count);
    void block_kick_one(int pos, int &alt, fp_t &p);

    virtual int alternate(int pos, fp_t fp) = 0; // get alternate position
    virtual int insert_to_bucket(fp_t *store, fp_t fp); // insert one fingerprint to bucket [pos] 
    virtual int lookup_in_bucket(int pos, fp_t fp); // lookup one fingerprint in  bucket [pos]
} ; 


template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::init(int _n, int _m, int _step)
{
    this -> n = _n;
    this -> m = _m;
    this -> max_kick_steps = _step;
    this -> filled_cell = 0;

    this -> block_number = ROUNDUP(int(_n / 0.95), 46) / 46;
    int how_many_bit = this -> block_number * 512;

    // here, we use memory_consumption as standard


    this -> memory_consumption = ROUNDUP(how_many_bit, 512) / 8; // how many bytes !
    this -> T = (uint8_t *) calloc(this -> memory_consumption, sizeof(char));

    this -> n = 64 * this -> block_number; // logical position number
    this -> m = 3; // m items perm bucket

    max_2_power = 1;

    for (; max_2_power * 2 < _n; ) max_2_power <<= 1;
}

template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::clear()
{
	this -> filled_cell = 0;
	//memset(this -> T, 0, sizeof(fp_t) * (this -> n * this -> m));
    memset(this -> T, 0, this -> memory_consumption);
}

template <typename fp_t, int fp_len>
fp_t MortonFilter<fp_t, fp_len>::fingerprint(int ele)
{
    fp_t h = HashUtil::BobHash32(&ele, 4) % ((1ull << fp_len) -1) + 1;
    return h;
}

template <typename fp_t, int fp_len>
int MortonFilter<fp_t, fp_len>::get_count(int block_number, int index)
{
    return (T[block_number * 512 / 8 + 368 / 8 + (index * 2) / 8] >> (index * 2 % 8)) & 3;
}

template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::set_count(int block_number, int index, int count)
{
    T[block_number * 512 / 8 + 368 / 8 + (index * 2) / 8] = ((uint8_t) T[block_number * 512 / 8 + 368 / 8 + (index * 2) / 8] & (uint8_t)(~(3 << (index * 2 % 8)))) + (count << (index * 2 % 8));
}

template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::print_block(int pos)
{
    int base_byte = pos / 64 * (512 / 8);
    int block_number = pos / 64;

    printf("block number = %d\n", block_number);
    printf("FCA : ");
    for (int i = 0; i < 64; i++)
        printf("%d ", get_count(block_number, i));
    puts("");

    printf("FSA : ");
    for (int i = 0; i < 46; i++)
        printf("%x ", T[base_byte + i]);
    puts("");
}

template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::get_bucket(int pos, fp_t *store)
{
    // T : uint8_t *
    int block = pos / 64;
    int index = pos % 64;
    int base_byte = block * 512 / 8;
    int offset_index = 0;

    for (int i = 0; i < index; i++) offset_index += get_count(block, i);

    /*
    if (pos == 1099) 
        print_block(pos);
    */

    int cur_bucket_count = get_count(block, index);

    for (int i = 0; i < 4; i++) store[i] = 0;
    for (int i = 0; i < cur_bucket_count; i++)
        store[i] = T[base_byte + offset_index + i];

    for (int i = index; i < 64; i++) offset_index += get_count(block, i);

    store[3] = offset_index; // store block load
}


template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::set_bucket(int pos, fp_t *store)
{
    int block = pos / 64;
    int index = pos % 64;
    int base_byte = block * 512 / 8;
    int offset_index = 0;

    for (int i = 0; i < index; i++) offset_index += get_count(block, i);

    int cur_bucket_count = get_count(block, index);
    int new_bucket_count = 0;

    for (int i = 0; i < 3; i++)
        if (store[i] != 0)
            new_bucket_count++;

    if (new_bucket_count < cur_bucket_count) // delete elements, move forward
    { 
        printf("delete pos = %d, index = %d, old = %d, new = %d\n", pos, index, cur_bucket_count, new_bucket_count);
        //print_block(pos);
        for (int i = offset_index + new_bucket_count; i + (cur_bucket_count - new_bucket_count) < 368 / 8; i++)
        {
            if (i + cur_bucket_count - new_bucket_count >= 368 / 8)
                T[base_byte + i] = 0;
            else
                T[base_byte + i] = T[base_byte + i + cur_bucket_count - new_bucket_count];
        }
    } else 
    if (new_bucket_count > cur_bucket_count) // add elements, move backward
    {
        //printf("add element to pos = %d, index = %d, old = %d, new = %d\n", pos, index, cur_bucket_count, new_bucket_count);
        for (int i = 368 / 8 - 1; i >= offset_index + new_bucket_count; i--)
            T[base_byte + i] = T[base_byte + i - new_bucket_count + cur_bucket_count];
    }

    set_count(block, index, new_bucket_count);
   // printf("%d 0x%x\n", base_byte + 368 / 8 + (index * 2) / 8, T[base_byte + 368 / 8 + (index * 2) / 8]);

    //printf("store : ");

    for (int i = 0; i < new_bucket_count; i++)
    {
        //printf("0x%x ", store[i]);
        T[base_byte + offset_index + i] = store[i];
    }
    //puts("");

    //print_block(base_byte);
}

/*
template <typename fp_t, int fp_len>
fp_t CuckooFilter<fp_t, fp_len>::get_item(int pos, int rk)
{
    return (fp_t) this -> T[pos * this -> m + rk];
}

template <typename fp_t, int fp_len>
void CuckooFilter<fp_t, fp_len>::set_item(int pos, int rk, fp_t fp)
{
    this -> T[pos * this -> m + rk] = (fp_t) fp;
}
*/
template <typename fp_t, int fp_len>
int MortonFilter<fp_t, fp_len>::insert_to_bucket(fp_t *store, fp_t fp)
{

    // if success return 0
    // if fail return 1

    //get_bucket(pos, store);

    if (store[3] >= 368 / 8) return 1; // block full !

    for (int i = 0; i < this -> m; i++) 
        if (store[i] == 0)
        {
            store[i] = fp;
            //set_bucket(pos, store);
            return 0;
        }

    return 1;
}

template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::block_kick_one(int pos, int &alt, fp_t &fp)
{
    //printf("ready to kick one\n");
    //print_block(pos);
    int block = pos / 64;
    int tmp = 0;
    for (int i = 63; i >= 0; i--)
    {
        if (get_count(block, i) > 0)
        {
            tmp = i;
            break;
        }
    }
    int old_count = get_count(block, tmp);
    set_count(block, tmp, old_count - 1);

    fp = T[block * 512 / 8 + 45];
    T[block * 512 / 8 + 45] = 0;
    alt = tmp + block * 64;

    //printf("\nkicked\n");
    //print_block(pos);
}
template <typename fp_t, int fp_len>
int MortonFilter<fp_t, fp_len>::insert(int ele)
{

    // If insert success return 0
    // If insert fail return 1

    fp_t fp = fingerprint(ele);

    int cur1 = this -> position_hash(ele);
    int cur2 = alternate(cur1, fp);
    if (alternate(cur2, fp) != cur1)
    {
        deln(ele);
        deln(fp);
        printf("cur1 = %d alt cur1 = %d\n", cur1, alternate(cur1, fp));
        printf("cur2 = %d alt cur2 = %d\n", cur2, alternate(cur2, fp));
    }

    fp_t store1[8];
    fp_t store2[8];

    /*
    if (ele == 0x6229393a)
    {
        print_block(cur1);
        print_block(cur2);
        fp = fingerprint(0x6229393a);
    }
    */
    int cnt1 = 0, cnt2 = 0;
    get_bucket(cur1, store1);
    for (int i = 0; i < this -> m; i++) cnt1 += (store1[i] != 0);
    get_bucket(cur2, store2);
    for (int i = 0; i < this -> m; i++) cnt2 += (store2[i] != 0);

    if (insert_to_bucket(store1, fp) == 0) {filled_cell++; set_bucket(cur1, store1); return 0;}
    if (insert_to_bucket(store2, fp) == 0) {filled_cell++; set_bucket(cur2, store2); return 0;}

    //get those item
    int cur;
    fp_t *cur_store;

    if (rand() & 1)
        cur = cur1, cur_store = store1;
    else 
        cur = cur2, cur_store = store2;

    fp_t tmp_fp = cur_store[0];
    if (tmp_fp == 0)
    {
        // this block is full ! kick another item in this block
        int kick_pos = 0;
        fp_t kick_item = 0;
        block_kick_one(cur, kick_pos, kick_item);

        cur_store[0] = fp;
        set_bucket(cur, cur_store);
        //printf("kick and store\n");
        //print_block(cur);

        cur = kick_pos;
        tmp_fp = kick_item;
    } else
    {
        cur_store[0] = fp;
        set_bucket(cur, cur_store);
    }

    int alt = alternate(cur, tmp_fp);
    
    for (int i = 0; i < this -> max_kick_steps; i++)
    {
        memset(store1, 0, sizeof(store1));
        get_bucket(alt, store1);
        if (insert_to_bucket(store1, tmp_fp) == 0) 
        {
            filled_cell++; 
            set_bucket(alt, store1);
            return 0;
        }

        fp = store1[0];
        if (fp == 0)
        {
            // this block is full ! kick another item in this block
            int kick_pos = 0;
            fp_t kick_item = 0;
            block_kick_one(alt, kick_pos, kick_item);

            store1[0] = tmp_fp;
            set_bucket(alt, store1);

            alt = kick_pos;
            fp = kick_item;
        } else
        {
            store1[0] = tmp_fp;
            set_bucket(alt, store1);
        }


        tmp_fp = fp;
        alt = alternate(alt, tmp_fp);
    }

    return 1;
}

template <typename fp_t, int fp_len>
int MortonFilter<fp_t, fp_len>::lookup_in_bucket(int pos, fp_t fp)
{
    // If lookup success return 1
    // If lookup fail and the bucket is full return 2
    // If lookup fail and the bucket is not full return 3

    fp_t store[8];
    get_bucket(pos, store);

    int isFull = 1;
    for (int i = 0; i < this -> m; i++)
    {
        fp_t t = store[i];
        if (t == fp) 
            return 1;
        isFull &= (t != 0);
    }
    return (isFull) ? 2 : 3;
}

template <typename fp_t, int fp_len>
bool MortonFilter<fp_t, fp_len>::lookup(int ele)
{

    // If ele is positive, return true
    // negative -- return false

	fp_t fp = fingerprint(ele);
    int pos1 = this -> position_hash(ele);
    int pos2 = alternate(pos1, fp);

    /*
    if (ele == 0x6229393a)
    {
        print_block(pos1);
        print_block(pos2);
        printf("fp = %x\n", fp);
    }
    */
  
    int ok1 = lookup_in_bucket(pos1, fp);

    if (ok1 == 1) return true;

    int ok2 = lookup_in_bucket(pos2, fp);

    return ok2 == 1;
}

template <typename fp_t, int fp_len>
double MortonFilter<fp_t, fp_len>::get_load_factor()
{
    return filled_cell * 1.0 / (block_number * 368 / 8);
}

template <typename fp_t, int fp_len>
void MortonFilter<fp_t, fp_len>::test_bucket()
{
    for (int i = 0; i < this -> n; i++)
    {
        fp_t store[8];
        store[0] = rand() % (1 << 16);
        store[1] = rand() % (1 << 16);
        store[2] = rand() % (1 << 16);

        set_bucket(i, store);

        fp_t tmp_store[8];
        get_bucket(i, tmp_store);
        for (int j = 0; j < 3; j++)
            assert(tmp_store[j] == store[j]);
    }
}

template <typename fp_t, int fp_len>
class StandardCuckooFilter : public SemiSortCuckooFilter<fp_t, fp_len>
{
    private : 
    int alternate(int pos, fp_t fp) // get alternate position
    {
        int ret =  pos ^ this->position_hash(HashUtil::MurmurHash32(fp));
        //if (ret == pos)
        //    puts("cuckoo xor 0");
        return ret;
    }
}; 

template <typename fp_t, int fp_len>
class XorFilter : public SemiSortCuckooFilter<fp_t, fp_len>
{
    private : 
    int alternate(int pos, fp_t fp) // get alternate position
    {
        int n = this -> n;

        int fp_hash = HashUtil::MurmurHash32(fp);

        //delta : the xor number
        int delta = this->position_hash(fp_hash);

        /*
        int ret;
        if (pos & 1)
            ret = pos - (delta | 1);
        else ret = pos + (delta | 1);
        return this -> position_hash(ret);
        */

        //if (delta == 0) delta = 1;

        //bias : the rotation bias
        int bias = this -> position_hash(HashUtil::MurmurHash32(fp_hash)) | 1;

        // add bias to avoid aggregation
        if (this -> debug_flag) printf("pos = %d, fp = %d, ", pos, fp);

        
        if (pos & 1)
            pos = pos + bias;
        else pos = pos - bias;
        
        pos = this -> position_hash(pos);
        //if (pos >= n) pos -= n;

        // find the corresponding segment of 'pos'
        // 1. pos ^ n 
        // ----- the highest different bit between position and n will be set to 1
        //
        // 2. find the highest bit 
        // ----- get the segment number
        //
        // 3. curlen = 1 << highest_bit
        // ----- get the segment length
        int segment_length = 1 << find_the_highest_bit(pos ^ n);

        // get the alternate position
        // 1. delta & (segment_length - 1)
        // ----- equals to delta % segment_length
        // 2. pos ^ ...
        // ----- xor (delta % segment_length)
        int t = (delta & (segment_length - 1));
        if (t == 0) t = 1;
        if (this -> debug_flag) printf("seg_len = %d, t = %d\n", segment_length, t);

        int ret = pos ^ t;

        // minus bias to avoid aggregation
        
        if (ret & 1)
            ret = ret + bias;
        else ret = ret - bias;
        
        ret = this -> position_hash(ret);

        return ret;
    }
}; 

template <typename fp_t, int fp_len>
class MortonAddFilter : public MortonFilter<fp_t, fp_len>
{
    private : 

    int alternate(int pos, fp_t fp)
    {
        // My cuckoo filter -- plus or minus
		//int n = this->n;
        //int B;
        //int bias = (n / (1 << fp_len)) * fp;  // bias -- avoid aggregate
		//pos = (pos + bias) % n;
        int delta = ((HashUtil::MurmurHash32(fp) & (this -> max_2_power - 1)) + 64) | 1;
        int ret = 0;

        if (pos & 1)
            ret = pos - delta;
        else ret = pos + delta;

        return this -> position_hash(ret);

    }
};

template <typename fp_t, int fp_len>
class RandomFilter : public SemiSortCuckooFilter<fp_t, fp_len>
{
    public : 
    
    int alt[1100][256];
    int id[1100];

    void init(int _n, int _m, int _steps)
    {
        SemiSortCuckooFilter<fp_t, fp_len>::init(_n, _m, _steps);
        for (int i = 0; i < _n; i++) id[i] = i;
        for (int i = 0; i < (1 << fp_len); i++)
        {
            random_shuffle(id, id + _n);
            for (int j = 0; j < _n; j += 2)
            {
                //printf("%d %d ", id[j], id[j + 1]);
                alt[id[j]][i] = id[j + 1], alt[id[j + 1]][i] = id[j];
            }
            //puts("");
        }
    }

    private : 

    int alternate(int pos, fp_t fp)
    {
        return alt[pos][fp];
    }
};

template <typename fp_t, int fp_len>
class AddSubFilter : public SemiSortCuckooFilter<fp_t, fp_len>
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
