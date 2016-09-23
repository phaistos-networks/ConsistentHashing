// In Consistent Hashing, a Ring  is represented as an array of sorted in ascending order tokens, and each of those
// tokens identifies a segment in the ring.
//
// The key property is that every segment owns the ring-space defined by the range:
// (prev_segment.token, segment.token]
// that is, starting but excluding the token of the previous segment in the ring, upto and including the token of the segment.
//
// The two exceptions are for tokens that <= the first tokens in the ring or > last tokens in the ring(ring semantics)
// -- For the last segment in the array, its next segment is the first segment in the array
// -- For the first segment in the array, its previous segment is the last segment in the array
//
//
// You should operate on ring segments, and for a typical distributed service, each segment will be owned by a primary replica, and based on
// replication strategies and factors, more(usually, the successor) segments will also get to hold to hold the same segment's data.
//
// You are probably going to be using a non-scalar type for the ring tokens (because you usually want at least 2^128 possible tokens, and so
// you want a nice struct that will represent those).
// To do so:
// You should provide template definitions for your own non-strandard token_t types like so:
// template<>
// static inline int8_t TrivialCmp<my_token_t>(const my_token_t &a, const my_token_t &b) { .. }
//
// and you should do the same for the minimum possible value for my_token_t like so:
// namespace std { template<> struct numeric_limits<my_token_t> { static const my_token_t min() { ... } }; }
// and you should likely need to do the same for std::min() and std::max() type specific implementations.
#pragma once
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <experimental/optional>

// Default implementation. Maybe a spaceship operator we could overload will make it into C++20?
template <typename T>
static inline int8_t TrivialCmp(const T &a, const T &b)
{
        return a == b ? 0 : (a < b ? -1 : 1);
}

namespace ConsistentHashing
{
	namespace
	{
		template<typename T>
		static void push_back_to(std::vector<T> *const out, const T*const values, const size_t cnt)
                {
                        for (uint32_t i{0}; i != cnt; ++i)
                                out->push_back(values[i]);
                }
        };

        // Returns lowest index where token <= tokens[index]
        // if it returns cnt, use 0 ( because tokens[0] owns ( tokens[cnt - 1], tokens[0] ]
        template <typename T>
        static uint32_t search(const T *const tokens, const uint32_t cnt, const T token)
        {
                int32_t h = cnt - 1, l{0};

                while (l <= h)
                {
                        // This protects from overflows: see http://locklessinc.com/articles/binary_search/
                        // The addition can be split up bitwise. The carry between bits can be obtained by
                        // the logical-and of the two summands. The resultant bit can be obtained by a XOR
                        //
                        // https://en.wikipedia.org/wiki/Binary_search_algorithm#Implementation_issues
                        // The problem with overflow is that if (l + h) add up to value greater than INT32_MAX,
                        // (exceeds the range of integers of the data type used to store the midpoint, even if
                        // l and h are withing rhe range). If l and h, this can non-negatives, this can be avoided
                        // by calculating the modpoint as: (l + (r - l) / 2)
                        //
                        // We are not using unsigned integers though -- though we should look for a way
                        // to do that so that we could safely use (l + (r - l ) / 2)
                        // so we can't use >> 1 here becuse (l + r) may result in a negative number
                        // and shifting by >> 1 won't divide that number by two.
                        const auto m = (l & h) + ((l ^ h) >> 1);
                        const auto v = tokens[m];
                        const auto r = TrivialCmp(token, v);

                        if (!r)
                                return m;
                        else if (r < 0)
                                h = m - 1;
                        else
                                l = m + 1;
                }

                return l;
        }

        // A segment in a ring. The segment is responsible(owns) the tokens range
        // (left, right] 	i.e left exlusive, right inclusive
        // whereas left is the token of the predecessor segment and right is the token of this segment
        template <typename token_t>
        struct ring_segment
        {
                token_t left;
                token_t right;

                ring_segment()
                {
                }

                ring_segment(const token_t l, const token_t r)
                    : left{l}, right{r}
                {
                }

                void set(const token_t l, const token_t r)
                {
                        left = l;
                        right = r;
                }

                // this segment's token
                auto token() const noexcept
                {
                        return right;
                }

                inline bool operator==(const ring_segment &o) const noexcept
                {
                        return left == o.left && right == o.right;
                }

                inline bool operator!=(const ring_segment &o) const noexcept
                {
                        return left != o.left || right != o.right;
                }

                inline bool operator<(const ring_segment &o) const noexcept
                {
                        return left < o.left || (left == o.left && right < o.right);
                }

                inline bool operator>(const ring_segment &o) const noexcept
                {
                        return left > o.left || (left == o.left && right > o.right);
                }

                int8_t cmp(const ring_segment &rhs) const noexcept
                {
                        if (tokens_wrap_around(left, right))
                        {
                                // there is only one segment that wraps around in the ring
                                return -1;
                        }
                        else if (tokens_wrap_around(rhs.left, rhs.right))
                        {
                                // there is only one segment that wraps around in the ring
                                return 1;
                        }
                        else
                        {
                                if (right == rhs.right)
                                        return 0;
                                else if (right > rhs.right)
                                        return 1;
                                else
                                        return -1;
                        }
                }

                [[gnu::always_inline]] static inline bool tokens_wrap_around(const token_t &l, const token_t &r) noexcept
                {
                        // true iff extends from last to the first ring segment
                        return l >= r;
                }

                bool contains(const ring_segment &that) const noexcept
                {
                        if (left == right)
                        {
                                // Full ring always contains all other ranges
                                return true;
                        }

                        const bool thisWraps = tokens_wrap_around(left, right);
                        const bool thatWraps = tokens_wrap_around(that.left, that.right);

                        if (thisWraps == thatWraps)
                                return left <= that.left && that.right <= right;
                        else if (thisWraps)
                        {
                                // wrapping might contain non-wrapping that is contained if both its tokens are in one of our wrap segments
                                return left <= that.left || that.right <= right;
                        }
                        else
                        {
                                // non-wrapping cannot contain wrapping
                                return false;
                        }
                }

                // For list of wrapped segments sorted by left token ascending, process the list to produce
                // an equivalent set of ranges, sans the overlapping ranges
                static void deoverlap(std::vector<ring_segment> *const segments)
                {
                        if (const auto cnt = segments->size())
                        {
                                const ring_segment *it = segments->data(), *const e = it + cnt;
                                ring_segment *out = segments->data(), cur = *it++;
                                const auto MinTokenValue = std::numeric_limits<token_t>::min();

                                while (it != e)
                                {
                                        // If current goes to the end of the ring, we are done
                                        if (cur.right == MinTokenValue)
                                        {
                                                // If one range is the full range, we return only that
                                                if (cur.left == MinTokenValue)
                                                {
                                                        segments->clear();
                                                        segments->push_back(cur);
                                                }
                                                else
                                                {
                                                        *out++ = ring_segment(cur.left, MinTokenValue);
                                                        segments->resize(out - segments->Values());
                                                }
                                                return;
                                        }

                                        const ring_segment next = *it++;

                                        // if the next left is == to the current right, we do not intersect per se, but replacing (A, B] and (B, C] by (A, C] is lefit,
                                        // and since this avoid special cases, and will result in more optimal segments, we transform here
                                        if (next.left <= cur.right)
                                        {
                                                // We do overlap
                                                if (next.right == MinTokenValue || cur.right < next.right)
                                                {
                                                        cur.right = next.right;
                                                }
                                        }
                                        else
                                        {
                                                *out++ = cur;
                                                cur = next;
                                        }
                                }
                                *out++ = cur;
                                segments->resize(out - segments->data());
                        }
                }

                // Copy of input list, with all segments unwrapped, sorted by left bound, and with overlapping bounds merged
                static void normalize(const ring_segment *const segments, const uint32_t segmentsCnt, std::vector<ring_segment> *const out)
                {
                        ring_segment res[2];

                        for (uint32_t i{0}; i != segmentsCnt; ++i)
                        {
                                if (const uint8_t n = segments[i].unwrap(res))
                                        push_back_to(out, res, n);
                        }

                        std::sort(out->begin(), out->end(), [](const auto &r1, const auto &r2) noexcept {
                                return r1.left < r2.left;
                        });

                        deoverlap(out);
                }

                static auto normalize(const ring_segment *const segments, const uint32_t segmentsCnt)
                {
                        std::vector<ring_segment> res;

                        normalize(segments, segmentsCnt, &res);
                        return res;
                }

                bool contains(const token_t &token) const noexcept
                {
                        if (wraps())
                        {
                                // We are wrapping around. Thee interval is (a, b] where a>= b
                                // then we have 3 cases which hold for any given token k, and we should return true
                                // 1. a < k
                                // 2. k <= b
                                // 3. b < k <= a
                                return token > left || right >= token;
                        }
                        else
                        {
                                // Range [a,b], a < b
                                return token > left && right >= token;
                        }
                }

                [[gnu::always_inline]] inline bool wraps() const noexcept
                {
                        return tokens_wrap_around(left, right);
                }

                inline bool intersects(const ring_segment that) const noexcept
                {
                        ring_segment out[2];

                        return intersection(that, out);
                }

                static uint8_t intersection_of_two_wrapping_segments(const ring_segment &first, const ring_segment &that, ring_segment *intersection) noexcept
                {
                        if (that.right > first.left)
                        {
                                intersection[0] = ring_segment(first.left, that.right);
                                intersection[1] = ring_segment(that.left, first.right);
                                return 2;
                        }
                        else
                        {
                                intersection[0] = ring_segment(that.left, first.right);
                                return 1;
                        }
                }

                static uint8_t intersection_of_single_wrapping_segment(const ring_segment &wrapping, const ring_segment &other, ring_segment *intersection) noexcept
                {
                        uint8_t size{0};

                        if (other.contains(wrapping.right))
                                intersection[size++] = ring_segment(other.left, wrapping.right);
                        if (other.contains(wrapping.left) && wrapping.left < other.right)
                                intersection[size++] = ring_segment(wrapping.left, other.right);

                        return size;
                }

                // Returns the intersection of two segments. That can be two disjoint ranges if one is wrapping and the other is not.
                // e.g for two nodes G and M, and a query range (D, T]; the intersection is (M-T] and (D-G]
                // If there is no interesection, an empty list is returned
                //
                // (12,7)^(5,20) => [(5,7), (12, 20)]
                // ring_segment(10, 100).intersection(50, 120) => [ ring_segment(50, 100) ]
                uint8_t intersection(const ring_segment &that, ring_segment *out) const noexcept
                {
                        if (that.contains(*this))
                        {
                                *out = *this;
                                return 1;
                        }
                        else if (contains(that))
                        {
                                *out = that;
                                return 1;
                        }
                        else
                        {
                                const bool thisWraps = tokens_wrap_around(left, right);
                                const bool thatWraps = tokens_wrap_around(that.left, that.right);

                                if (!thisWraps && !thatWraps)
                                {
                                        // Neither wraps; fast path
                                        if (!(left < that.right && that.left < right))
                                                return 0;

                                        *out = ring_segment(std::max<token_t>(left, that.left), std::min<token_t>(right, that.right));
                                        return 1;
                                }
                                else if (thisWraps && thatWraps)
                                {
                                        // Two wrapping ranges always intersect. 
					// We have determined that neither this or that contains the other, we are left
                                        // with two possibilities and mirror images of each such case:
                                        // 1. both of s (1,2] endpoints lie in this's (A, B] right semgnet:
                                        // 2. only that's start endpoint lies in this's right segment:
                                        if (left < that.left)
                                                return intersection_of_two_wrapping_segments(*this, that, out);
                                        else
                                                return intersection_of_two_wrapping_segments(that, *this, out);
                                }
                                else if (thisWraps && !thatWraps)
                                        return intersection_of_single_wrapping_segment(*this, that, out);
                                else
                                        return intersection_of_single_wrapping_segment(that, *this, out);
                        }
                }

                // Subtracts a portion of this range
                // @contained : The range to subtract from `this`: must be totally contained by this range
                // @out: List of ranges left after subtracting contained from `this` (@return value is size of @out)
                //
                // i.e ring_segment(10, 100).subdvide(ring_segment(50, 55)) => [ ring_segment(10, 50), ring_segment(55, 110) ]
                uint8_t subdivide(const ring_segment &contained, ring_segment *const out) const noexcept
                {
                        uint8_t size{0};

                        if (left != contained.left)
                                out[size++] = ring_segment(left, contained.left);
                        if (right != contained.right)
                                out[size++] = ring_segment(contained.right, right);
                        return size;
                }

		// if this segment wraps, it will return two segments
		// 1. (left, std::numeric_limits<token_t>::min())
		// 2. (std::numeric_limits<token_t>::min(), right)
		// otherwise, it will return itself
                uint8_t unwrap(ring_segment *const out) const noexcept
                {
                        const auto MinTokenValue = std::numeric_limits<token_t>::min();

                        if (false == wraps() || right == MinTokenValue)
                        {
                                *out = *this;
                                return 1;
                        }
                        else
                        {
                                out[0] = ring_segment(left, MinTokenValue);
                                out[1] = ring_segment(MinTokenValue, right);
                                return 2;
                        }
                }

                // Compute difference betweet two ring segments
                // This is very handy for computing, e.g the segments a node will need to fetch, when moving to a given token
                // e.g segment(5, 20).difference(segment(2, 25)) => [ (2, 5), (20, 25) ]
                // e.g segment(18, 25).difference(segment(5,20)) => [ (5, 18) ]
                //
                // In other words, compute the missing segments(ranges) that (*this) is missing from rhs
                uint8_t difference(const ring_segment &rhs, ring_segment *const result) const
                {
                        ring_segment intersectionSet[2];

                        switch (intersection(rhs, intersectionSet))
                        {
                                case 0:
                                        // not intersected
                                        *result = rhs;
                                        return 1;

                                case 1:
                                        // compute missing sub-segments
                                        return rhs.subdivide(intersectionSet[0], result);

                                default:
                                {
                                        const auto first = intersectionSet[0], second = intersectionSet[1];
                                        ring_segment tmp[2];

                                        rhs.subdivide(first, tmp);
                                        // two intersections; subtracting only one of them will yield a single segment
                                        return tmp[0].subdivide(second, result);
                                }
                        }
                }

                // split the segment in two, halved at segmentToken value (if segmentToken it makes sense)
                //
                // i.e ring_segment(10, 20).split(18) => (  ring_segment(10, 18), ring_segment(18, 20) )
                std::experimental::optional<std::pair<ring_segment, ring_segment>> split(const token_t segmentToken) const noexcept
                {
                        if (left == segmentToken || right == segmentToken || !contains(segmentToken))
                                return {};

                        return {{ring_segment(left, segmentToken), ring_segment(segmentToken, right)}};
                }


		// Two handy utility methods for checking if a token is in 0+ segments
                static bool token_in_any_segments(const token_t token, const std::vector<ring_segment> &segments)
                {
                        const auto *const r = segments.data();
                        const auto n = segments.size();

                        switch (n)
                        {
                                case 4:
                                        if (r[3].contains(token))
                                                return true;
                                case 3:
                                        if (r[2].contains(token))
                                                return true;
                                case 2:
                                        if (r[1].contains(token))
                                                return true;
                                case 1:
                                        if (r[0].contains(token))
                                                return true;
                                case 0:
                                        break;

                                default:
                                        for (uint32_t i{0}; i != n; ++i)
                                        {
                                                if (r[i].contains(token))
                                                        return true;
                                        }
                                        break;
                        }

                        return false;
                }

		// if input segments are ordered, you should consider using this method instead
                static bool token_in_any_ordered_segments_with_base(const token_t token, const std::vector<ring_segment> &segments, uint32_t &base)
                {
                        const auto *const all = segments.data();
                        const auto cnt = segments.size();

                        if (base == 0 && token < all[0].left && all[cnt - 1].wraps())
                                return token <= all[cnt - 1].right;
                        else
                        {
                                for (uint32_t i = base; i < cnt; ++i)
                                {
                                        if (token > all[i].right)
                                                ++base;
                                        else
                                                return all[i].contains(token);
                                }

                                return false;
                        }
                }
        };


	// A Ring of tokens
        template <typename T>
        struct Ring
        {
                using token_t = T;
		using segment_t = ring_segment<T>;

                const T *const tokens;
                const uint32_t cnt;

                Ring(const T *const v, const uint32_t n)
                    : tokens{v}, cnt{n}
                {
                }

                inline auto size() const noexcept
                {
                        return cnt;
                }

                uint32_t index_of(const T token) const noexcept
                {
                        for (int32_t h = cnt - 1, l{0}; l <= h;)
                        {
                                const auto m = (l & h) + ((l ^ h) >> 1);
                                const auto v = tokens[m];
                                const auto r = TrivialCmp(token, v);

                                if (!r)
                                        return m;
                                else if (r < 0)
                                        h = m - 1;
                                else
                                        l = m + 1;
                        }

                        return UINT32_MAX;
                }

                inline bool is_set(const T token) const noexcept
                {
                        return index_of(token) != UINT32_MAX;
                }

                inline uint32_t search(const T token) const noexcept
                {
                        return ConsistentHashing::search(tokens, cnt, token);
                }

                // In a distributed systems, you should map the token to a node (or the segment index returned by this method)
                inline uint32_t index_owner_of(const T token) const noexcept
                {
                        // modulo is not cheap, and comparisons are much cheaper, but branchless is nice
                        return search(token) % cnt;
                }

                inline auto token_owner_of(const T token) const noexcept
                {
                        return tokens[index_owner_of(token)];
                }

                const T &token_predecessor(const T token) const noexcept
                {
                        const auto idx = index_of(token);

                        return tokens[(idx + (cnt - 1)) % cnt];
                }

                const T &token_successor(const T token) const noexcept
                {
                        const auto idx = index_of(token);

                        return tokens[(idx + 1) % cnt];
                }

                // segment in the ring that owns this token
                // based on the (prev segment.token, this segment.token] ownership rule
                auto token_segment(const T token) const noexcept
                {
                        const auto idx = index_of(token);

                        return ring_segment<T>(tokens[(idx + (cnt - 1)) % cnt], tokens[idx]);
                }

                void segments(std::vector<ring_segment<T>> *const res) const
                {
                        if (cnt)
                        {
                                res->reserve(cnt + 2);
                                for (uint32_t i{1}; i != cnt; ++i)
                                        res->push_back({tokens[i - 1], tokens[i]});
                                res->push_back({tokens[cnt - 1], tokens[0]});
                        }
                }

                auto segments() const
                {
                        std::vector<ring_segment<T>> res;

                        segments(&res);
                        return res;
                }

		// Assuming a node is a replica for tokens in segments `current`, and then the segments is is a replica for
		// change to `updated`. 
		//
		// This handy utility method will generate a pair of segments list. 
		// The first is segments the node will need to *fetch* from other nodes in the ring, because it will now be also responsible
		// for those segments, but it does not have the data.
		// The second is segments the node will need to *stream* to other nodes in the ring, because it will no longer hold data for them.
		//
		// Obviously, if a node is just introduced to a ring (i.e have only updated and no current segments ), it should
		// just fetch data for all the current segments. Conversely, if the node is exiting the ring, it should
		// consider streaming all the data it has to other nodes if needed, and not fetch any data to itself.
		static auto compute_segments_ownership_updates(const std::vector<segment_t> &current, const std::vector<segment_t> &updated)
                {
			std::vector<segment_t> toFetch, toStream;
                        segment_t segmentsList[2];

			for (const auto curSegment : current)
                        {
                                bool intersect{false};

				for (const auto updatedSegment : updated)
                                {
                                        if (curSegment.intersects(updatedSegment))
                                        {
                                                push_back_to(&toStream, segmentsList, updatedSegment.difference(curSegment, segmentsList));
                                                intersect = true;
                                        }
                                }

                                if (!intersect)
                                        toStream.push_back(curSegment); 
                        }

			for (const auto updatedSegment : updated)
                        {
                                bool intersect{false};

				for (const auto curSegment : current)
                                {
                                        if (updatedSegment.intersects(curSegment))
                                        {
                                                push_back_to(&toFetch, segmentsList, curSegment.difference(updatedSegment, segmentsList));
                                                intersect = true;
                                        }
                                }

                                if (!intersect)
                                        toFetch.push_back(updatedSegment);
                        }

			return std::make_pair(toFetch, toStream);
                }
        };
}
