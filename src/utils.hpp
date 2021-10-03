#ifndef UTILS_RING
#define UTILS_RING


// A hashing function on pairs of uint64_t

namespace std {
    template <> struct hash<std::pair<uint64_t, uint64_t>> {
        inline size_t operator()(const std::pair<uint64_t, uint64_t> &v) const {
            std::hash<uint64_t> int_hasher;
            return int_hasher(v.first) ^ int_hasher(v.second);
        }
    };
}

uint32_t* dat;

// Utility function for suffix array construction
int compare (const void *p1, const void *p2)
{
    uint64_t r1 = *((uint64_t*)p1);
    uint64_t r2 = *((uint64_t*)p2);
    if (r1 == r2) return 0;
    while (dat[r1] == dat[r2]) { r1++; r2++; }
    return dat[r1]-dat[r2];
}

#endif
