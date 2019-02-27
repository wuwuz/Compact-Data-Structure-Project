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

    void init(int _n, int _m, int _max_kick_steps)
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
} ; 


template <typename fp_t, int fp_len>
void CuckooFilter<fp_t, fp_len>::init(int _n, int _m, int _step)
{
    this -> n = _n;
    this -> m = _m;
    this -> max_kick_steps = _step;
    this -> filled_cell = 0;
    this -> memory_consumption = _n * _m * sizeof(fp_t);

    max_2_power = 1;
    for (; max_2_power * 2 < _n; ) max_2_power <<= 1;

    this -> T = (fp_t*) calloc(this -> n * this -> m, sizeof(fp_t));
}

template <typename fp_t, int fp_len>
void CuckooFilter<fp_t, fp_len>::clear()
{
	this -> filled_cell = 0;
	memset(this -> T, 0, sizeof(fp_t) * (this -> n * this -> m));
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

    uint32_t *T;
    uint32_t encode_table[1 << 16];
    uint32_t decode_table[1 << 16];

	~SemiSortCuckooFilter() { free(T); }
    
    int filled_cell;
    int max_kick_steps;

    fp_t fingerprint(int ele); // 32-bit to 'fp_len'-bit fingerprint

    //interface for semi-sorted bucket
    void get_bucket(int pos, fp_t *store);
    void set_bucket(int pos, fp_t *sotre);
    void test_bucket();

    virtual int alternate(int pos, fp_t fp) = 0; // get alternate position
    virtual int insert_to_bucket(fp_t *store, fp_t fp); // insert one fingerprint to bucket [pos] 
    virtual int lookup_in_bucket(int pos, fp_t fp); // lookup one fingerprint in  bucket [pos]
} ; 


template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::init(int _n, int _m, int _step)
{
    this -> n = _n;
    this -> m = _m;
    this -> max_kick_steps = _step;
    this -> filled_cell = 0;

    int how_many_bit = this -> n * this -> m * (fp_len - 1);

    this -> memory_consumption = ROUNDUP(ROUNDUP(how_many_bit, 8) / 8 + 8, 8); // how many bytes !

    max_2_power = 1;
    for (; max_2_power * 2 < _n; ) max_2_power <<= 1;

    //this -> T = (fp_t*) calloc(this -> n * this -> m, sizeof(fp_t));
    this -> T = (uint32_t *) calloc(this -> memory_consumption, sizeof(char));

    int index = 0;
    for (int i = 0; i < 16; i++)
        for (int j = 0; j <= i; j++)
            for (int k = 0; k <= j; k++)
                for (int l = 0; l <= k; l++)
                {
                    int plain_bit = (i << 12) + (j << 8) + (k << 4) + l;
                    encode_table[plain_bit] = index;
                    decode_table[index] = plain_bit;
                    ++index;
                }
    deln(index);
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
int SemiSortCuckooFilter<fp_t, fp_len>::insert_to_bucket(fp_t *store, fp_t fp)
{

    // if success return 0
    // if fail return 1

    //get_bucket(pos, store);

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
    if (ok1 == 3) return false;

    int ok2 = lookup_in_bucket(pos2, fp);

    return ok2 == 1;
}

template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_load_factor()
{
    return filled_cell * 1.0 / this -> n / this -> m;
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::test_bucket()
{
    for (int i = 0; i < this -> n; i++)
    {
        fp_t store[8];
        store[0] = rand() % (1 << 16);
        store[1] = rand() % (1 << 16);
        store[2] = rand() % (1 << 16);
        store[3] = rand() % (1 << 16);

        set_bucket(i, store);

        fp_t tmp_store[8];
        get_bucket(i, tmp_store);
        for (int j = 0; j < 4; j++)
            assert(tmp_store[j] == store[j]);
    }
}

template <typename fp_t, int fp_len>
class StandardCuckooFilter : public CuckooFilter<fp_t, fp_len>
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
        if (delta == 0) delta = 1;

        //bias : the rotation bias
        int bias = delta;

        // add bias to avoid aggregation
        pos = pos + bias;
        if (pos >= n) pos -= n;

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
        int alternate_position = pos ^ (delta & (segment_length - 1));

        // minus bias to avoid aggregation
        alternate_position = alternate_position - bias;
        if (alternate_position < 0) alternate_position += n;

        return alternate_position;
    }
}; 

template <typename fp_t, int fp_len>
class AddFilter : public CuckooFilter<fp_t, fp_len>
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
class MyFilter : public CuckooFilter<fp_t, fp_len>
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
