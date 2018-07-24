#include <bits/stdc++.h>
#include <cstdint>
#define memcle(a) memset(a, 0, sizeof(a))
#define sqr(a) ((a) * (a))
#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a) cerr << #a << " = " << a << endl
typedef long long LL;
using namespace std;
const int N = 2000100;
const int K = 200; // threshold for block split

int exam[N];

/**
 * @brief mix 3 32-bit values reversibly
 *
 * For every delta with one or two bits set, and the deltas of all three
 * high bits or all three low bits, whether the original value of a,b,c
 * is almost all zero or is uniformly distributed.
 *
 * If mix() is run forward or backward, at least 32 bits in a,b,c
 * have at least 1/4 probability of changing.
 *
 * If mix() is run forward, every bit of c will change between 1/3 and
 * 2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
 * mix() was built out of 36 single-cycle latency instructions in a 
 * structure that could supported 2x parallelism, like so:
 *    a -= b; 
 *    a -= c; x = (c>>13);
 *    b -= c; a ^= x;
 *    b -= a; x = (a<<8);
 *    c -= a; b ^= x;
 *    c -= b; x = (b>>13);
 *     ...
 *
 * Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
 * of that parallelism.  They've also turned some of those single-cycle
 * latency instructions into multi-cycle latency instructions. Still,
 * this is the fastest good hash I could find. There were about 2^68
 * to choose from. I only looked at a billion or so.
 */
#define BOBHASH_MIX(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8);  \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12); \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5);  \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/**
 * Every bit of the key affects every bit of the return value. 
 * Every 1-bit and 2-bit delta achieves avalanche.
 * About 6 * length + 35 instructions.
 *
 * Use for hash table lookup, or anything where one collision in 2^32 is acceptable.
 * Do NOT use for cryptographic purposes.
 */
  uint32_t
bobhash (const void *key, size_t key_size)
{
  const uint32_t BOBHASH_GOLDEN_RATIO = 0x9e3779b9;
  uint32_t a = BOBHASH_GOLDEN_RATIO;
  uint32_t b = BOBHASH_GOLDEN_RATIO;
  uint32_t c = 0;
  uint32_t length = key_size;

  uint8_t* work_key = (uint8_t*) key;

  /* handle most of the key */
  while (length >= 12)
  {
    a += (work_key[0] + ((uint32_t)work_key[1] << 8) + ((uint32_t)work_key[2] << 16) + ((uint32_t)work_key[3] << 24));
    b += (work_key[4] + ((uint32_t)work_key[5] << 8) + ((uint32_t)work_key[6] << 16) + ((uint32_t)work_key[7] << 24));
    c += (work_key[8] + ((uint32_t)work_key[9] << 8) + ((uint32_t)work_key[10] << 16)+ ((uint32_t)work_key[11] << 24));
    BOBHASH_MIX (a,b,c);
    work_key += 12; 
    length -= 12;
  }

  /* handle the last 11 bytes */
  c += key_size;
  switch (length)
  {
    case 11: c += ((uint32_t)work_key[10] << 24);
    case 10: c += ((uint32_t)work_key[9] << 16);
    case 9 : c += ((uint32_t)work_key[8] << 8);
    case 8 : b += ((uint32_t)work_key[7] << 24);
    case 7 : b += ((uint32_t)work_key[6] << 16);
    case 6 : b += ((uint32_t)work_key[5] << 8);
    case 5 : b += work_key[4];
    case 4 : a += ((uint32_t)work_key[3] << 24);
    case 3 : a += ((uint32_t)work_key[2] << 16);
    case 2 : a += ((uint32_t)work_key[1] << 8);
    case 1 : a += work_key[0];
  }

  BOBHASH_MIX (a,b,c);

  return c & (0x7fffffff); // 32-bit
}

// 2-4 Cuckoo Filter
struct Filter
{
    int n; // number of buckets
    int m; // number of slots per bucket
    int filled_cell;
    char *T; // table
    int max_kick_steps;

    int position_hash(int ele)
    {
        return ele % n;
    }

    int Hash_32(int h)
    {
        h ^= h >> 16;
        h *= 0x85ebca6b;
        h ^= h >> 13;
        h *= 0xc2b2ae35;
        h ^= h >> 16;
        return h;
    }

    int initialization(int _n, int _m)
    {
        n = _n;
        m = _m;
        max_kick_steps = 200;
        filled_cell = 0;

        T = (char*) calloc(n * m, sizeof(char));
    } 

    char fingerprint(int ele)
    {
        char h = bobhash(&ele, 8) % 255 + 1;
        return h;
    }

    char get_item(int pos, int rk)
    {
        return (char) T[pos * m + rk];
    }

    int set_item(int pos, int rk, char fp)
    {
        T[pos * m + rk] = (char) fp;
    }

    int alternate(int pos, char fp)
    {
        // standard cuckoo filter -- xor
        return pos ^ position_hash(Hash_32(fp));
    }

    int insert_to_bucket(char fp, int pos)
    {

        // if success return 0
        // if fail return 1

        //debug(fp);
        //deln(pos);

        for (int i = 0; i < m; i++) 
            if (get_item(pos, i) == 0)
            {
                set_item(pos, i, fp);
                return 0;
            } else 
            {
                char t = get_item(pos, i);
                //deln(t);
            }

        return 1;
    }

    int insert(int ele)
    {

        // If insert success return 0
        // If insert fail return 1
        //
        //debug(ele);

        char fp = fingerprint(ele);
        int cur1 = position_hash(ele);
        int cur2 = alternate(cur1, fp);

        //debug(fp);
        //debug(cur1);
        //deln(cur2);

        for (int i = 0; i < max_kick_steps; i++)
            if (insert_to_bucket(fp, cur1) == 1) 
            {
                if (insert_to_bucket(fp, cur2) == 0) 
                {
                    filled_cell++;
                    return 0;
                }

                int cur = cur1;
                if (rand() % 2 == 0) cur = cur2;

                // fail -- this bucket is full
                int rk = rand() % m;
                char tmp_fp = get_item(cur, rk);
                set_item(cur, rk, fp);

                fp = tmp_fp;
                cur1 = cur;
                cur2 = alternate(cur1, fp);
            } else 
            {
                filled_cell++;
                return 0; // insert success!
            }

        return 1;
    }

    bool lookup_in_bucket(int pos, char fp)
    {
        // If lookup success return true
        // If lookup fail return false

        for (int i = 0; i < m; i++)
        {
            char t = get_item(pos, i);
            if (get_item(pos, i) == (char) fp) 
                return true;
        }
        return false;
    }

    bool lookup(int ele)
    {

        // If ele is positive, return true
        // negative -- return false

        char fp = fingerprint(ele);
        int pos1 = position_hash(ele);
        int pos2 = alternate(pos1, fp);

        //debug(ele); debug(fp); debug(pos1); debug(pos2);

        bool ok1 = lookup_in_bucket(pos1, fp);
        bool ok2 = lookup_in_bucket(pos2, fp);

        //debug(ok1); deln(ok2);

        return ok1 | ok2;
    }

} a;

int main()
{

    srand(12345);

    int max_item = 1000000;

    int n = 1;
    //printf("%d %d\n", bobhash(&n, 32), bobhash(&n, 32));
    for (; n * 4 <= max_item; n <<= 1);
    int m = 4;

    a.initialization(n, m);

    for (int i = 1; i <= max_item * 2; i++) exam[i] = i;
    random_shuffle(exam + 1, exam + 1 + max_item * 2);

    int fail_insert = 0;
    int false_positive = 0;
    int false_negative = 0;

    for (int i = 1; i <= max_item; i++) 
    {
        if (a.insert(exam[i]) == 1) // insert fail
            fail_insert++;
    }

    // positive lookup

    for (int i = 1; i <= max_item; i++)
    {
        if (a.lookup(exam[i]) == false) // false neage
        {
            //deln(i);
            false_negative++;
        }
    }

    for (int i = max_item + 1; i <= max_item * 2; i++)
    {
        if (a.lookup(exam[i]) == true) // false positive
            false_positive++;
    }

    printf("testing : n = %d, m = %d\n", n, m);
    printf("insert item = %d\n", max_item);
    printf("insert fail rate = %.2f\n", fail_insert * 100.0 / max_item);
    printf("false negative = %.2f\n", false_negative * 100.0 / max_item);
    printf("false positive = %.2f\n", false_positive * 100.0 / max_item);
    printf("load factor = %.2f\n", a.filled_cell * 100.0 / m / n);

    return 0;
}

