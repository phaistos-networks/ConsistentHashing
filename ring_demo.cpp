#include "consistent_hashing.h"
#include <stdio.h>

int main()
{
#if 0
        uint32_t ringTokens[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
        const ConsistentHashing::Ring<uint32_t> ring(ringTokens, sizeof(ringTokens) / sizeof(ringTokens[0]));

        printf("%u\n", ring.search(12));
#endif

	using segment = ConsistentHashing::ring_segment<uint32_t>;
        std::vector<segment> current({{5, 10}, {11, 20}, {21, 30}, {31, 40}, {41, 50}});
        std::vector<segment> updated({{2, 8}, {45, 110}});
        const auto res = ConsistentHashing::Ring<uint32_t>::compute_segments_ownership_updates(current, updated);

	printf("%u %u\n", res.first.size(), res.second.size());
	return 0;

}
