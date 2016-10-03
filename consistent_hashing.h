// In Consistent Hashing, a Ring is represented as an array of sorted in ascending order tokens, and each of those
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
#pragma once
#ifdef HAVE_SWITCH
#include "switch.h"
#include "switch_vector.h"
#include <experimental/optional>
#else
#include <algorithm>
#include <experimental/optional>
#include <limits>
#include <stdint.h>
#include <string.h>
#include <vector>
#endif
#include <unordered_map>

namespace ConsistentHashing
{
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
        // See also: https://en.wikipedia.org/wiki/Circular_segment
        template <typename token_t>
        struct ring_segment
        {
                token_t left;
                token_t right;

                uint64_t span() const noexcept
                {
                        if (wraps())
                        {
                                require(left >= right);
                                return uint64_t(std::numeric_limits<token_t>::max()) - left + right;
                        }
                        else
                        {
                                require(right >= left);
                                return right - left;
                        }
                }

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

                // masks a segment `mask` from a segment `s`, if they intersect, and return 0+ segments
                //
                // It is very important that we get this right, otherwise other methods that depend on it will produce crap
                // returns a pair, where the first is true if the segment was intersected by the mask, false otherwise, and the second
                // is the number of segments it was partitioned to (can be 0)
                std::pair<bool, uint8_t> mask(const ring_segment mask, ring_segment *const out) const noexcept
                {
                        if (false == intersects(mask))
                                return {false, 0};
                        else if (mask.contains(*this))
                        {
                                // completely masked
                                return {true, 0};
                        }
                        else
                        {
                                // partially masked
                                uint8_t n{0};

                                if (mask.wraps() || wraps())
                                        n = mask.difference(*this, out);
                                else if (mask.right > left)
                                {
                                        if (mask.left < right && mask.left > left)
                                                out[n++] = {left, mask.left};

                                        if (mask.right < right)
                                                out[n++] = {mask.right, right};
                                }

                                return {true, n};
                        }
                }

                static void mask_segments_impl(const ring_segment *it, const ring_segment *const end, const std::vector<ring_segment> &toExclude, std::vector<ring_segment> *const out)
                {
                        ring_segment list[2];

                        for (auto i{it}; i != end; ++i)
                        {
                                const auto in = *i;

                                for (const auto mask : toExclude)
                                {
                                        if (const auto res = in.mask(mask, list); res.first)
                                        {
                                                // OK, either completely or partially masked

                                                if (res.second)
                                                        mask_segments_impl(list, list + res.second, toExclude, out);

                                                goto next;
                                        }
                                }

                                out->push_back(in);

                        next:;
                        }
                }

                static void mask_segments(const ring_segment *it, const ring_segment *const end, const std::vector<ring_segment> &toExclude, std::vector<ring_segment> *const out)
                {
                        if (toExclude.size())
                        {
                                mask_segments_impl(it, end, toExclude, out);
                                // Just in case (this is cheap)
                                sort_and_deoverlap(out);
                        }
                        else
                                out->insert(out->end(), it, end);
                }

                static void mask_segments(const std::vector<ring_segment> &in, const std::vector<ring_segment> &toExclude, std::vector<ring_segment> *const out)
                {
                        mask_segments(in.data(), in.data() + in.size(), toExclude, out);
                }

                static auto mask_segments(const std::vector<ring_segment> &in, const std::vector<ring_segment> &toExclude)
                {
                        std::vector<ring_segment> out;

                        mask_segments(in.begin(), in.end(), &out);
                        return out;
                }

                // For list of wrapped segments sorted by left token ascending, process the list to produce
                // an equivalent set of ranges, sans the overlapping ranges
                // it will also merge together ranges
                // i.e [(8, 10],(8, 15],(14, 18],(17, 18]] => [ (8, 18] ]
                //
                // this will only work if the segments are properly sorted. see sort_and_deoverlap()
                // This utility method deals with invalid segments as well (e.g you can't really have more than one segments that wrap)
                static void deoverlap(std::vector<ring_segment> *const segments)
                {
                        auto out = segments->data();

                        for (auto *it = segments->data(), *const end = it + segments->size(); it != end;)
                        {
                                auto s = *it;

                                if (it->right <= it->left)
                                {
                                        // This segment wraps
                                        // deal with e.g [30, 4], [35, 8], [40, 2]
                                        // that'd be an invalid list of segments(there can only be one wrapping segment), but we 'll deal with it anyway
                                        const auto wrappedSegmentIt = it;

                                        for (++it; it != end; ++it)
                                        {
                                                if (it->right > s.right)
                                                        s.right = it->right;
                                        }

                                        // we need to potentially drop some of them segments if the wrapping segment overlaps them
                                        if (wrappedSegmentIt != (it = segments->data()) && s.right >= it->right)
                                        {
                                                s.right = it->right;
                                                memmove(it, it + 1, (out - it) * sizeof(ring_segment));
                                                --out;
                                        }

                                        *out++ = s;
                                        break;
                                }
                                else
                                {
                                        for (++it; it != end && ((*it == s) || (it->left >= s.left && s.right > it->left)); ++it)
                                                s.right = it->right;

                                        if (out == segments->data() || false == out[-1].contains(s))
                                        {
                                                // deal with (8, 30],(9, 18]
                                                *out++ = s;
                                        }
                                }
                        }

                        segments->resize(out - segments->data());

                        if (segments->size() == 1 && segments->back().left == segments->back().right)
                        {
                                // spans the whole ring
                                const auto MinTokenValue = std::numeric_limits<token_t>::min();

                                segments->pop_back();
                                segments->push_back({MinTokenValue, MinTokenValue});
                        }
                }

                // utility method; sorts segments so that deoverlap() can process them
                static void sort_and_deoverlap(std::vector<ring_segment> *const segments)
                {
                        std::sort(segments->begin(), segments->end(), [](const auto &a, const auto &b) { return a.left < b.left || (a.left == b.left && a.right < b.right); });
                        deoverlap(segments);
                }

                // Copy of input list, with all segments unwrapped, sorted by left bound, and with overlapping bounds merged
                static void normalize(const ring_segment *const segments, const uint32_t segmentsCnt, std::vector<ring_segment> *const out)
                {
                        ring_segment res[2];

                        for (uint32_t i{0}; i != segmentsCnt; ++i)
                        {
                                if (const uint8_t n = segments[i].unwrap(res))
                                        out->insert(out->end(), res, res + n);
                        }

                        sort_and_deoverlap(out);
                }

                static auto normalize(const ring_segment *const segments, const uint32_t segmentsCnt)
                {
                        std::vector<ring_segment> res;

                        normalize(segments, segmentsCnt, &res);
                        return res;
                }

                // true iff segment contains the token
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

                static uint8_t _intersection_of_two_wrapping_segments(const ring_segment &first, const ring_segment &that, ring_segment *intersection) noexcept
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

                static uint8_t _intersection_of_single_wrapping_segment(const ring_segment &wrapping, const ring_segment &other, ring_segment *intersection) noexcept
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
                // see also mask()
                //
                // this is the result of the logical operation: ((*this) & that)
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
                                        // 1. both of s (1,2] endpoints lie in this's (A, B] right segment
                                        // 2. only that's start endpoint lies in this's right segment:
                                        if (left < that.left)
                                                return _intersection_of_two_wrapping_segments(*this, that, out);
                                        else
                                                return _intersection_of_two_wrapping_segments(that, *this, out);
                                }
                                else if (thisWraps && !thatWraps)
                                        return _intersection_of_single_wrapping_segment(*this, that, out);
                                else
                                        return _intersection_of_single_wrapping_segment(that, *this, out);
                        }
                }

                // Subtracts a portion of this range
                // @contained : The range to subtract from `this`: must be totally contained by this range
                // @out: List of ranges left after subtracting contained from `this` (@return value is size of @out)
                //
                // i.e ring_segment(10, 100).subdvide(ring_segment(50, 55)) => [ ring_segment(10, 50), ring_segment(55, 110) ]
                //
                // You may want to use mask() instead, which is more powerful and covers wrapping cases, etc
                uint8_t subdivide(const ring_segment &contained, ring_segment *const out) const noexcept
                {
                        if (contained.contains(*this))
                        {
                                // contained actually contains this segment
                                return 0;
                        }

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
                // There is an opposite operation, mask()
                //
                // This is the result of the logical operation: (rhs & (~(rhs & (*this))) )
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

                // split the segment in two, halved at segmentToken value (if segmentToken is contained in segment)
                //
                // i.e ring_segment(10, 20).split(18) => (  ring_segment(10, 18), ring_segment(18, 20) )
                std::experimental::optional<std::pair<ring_segment, ring_segment>> split(const token_t segmentToken) const noexcept
                {
                        if (left == segmentToken || right == segmentToken || !contains(segmentToken))
                                return {};

                        return {{ring_segment(left, segmentToken), ring_segment(segmentToken, right)}};
                }

#ifdef HAVE_SWITCH
                void serialize(IOBuffer *const b) const
                {
                        b->Serialize(left);
                        b->Serialize(right);
                }

                void deserialize(ISerializer *const b) const
                {
                        b->Unserialize<token_t>(&left);
                        b->Unserialize<token_t>(&right);
                }
#endif

                // Make sure segments is properly ordered and deoverlapped
                // see sort_and_deoverlap()
                static bool segments_contain(const token_t token, const ring_segment *const segments, const uint32_t cnt)
                {
                        if (!cnt)
                                return false;

                        int32_t h = cnt - 1;

                        if (segments[h].wraps())
                        {
                                if (segments[h--].contains(token))
                                {
                                        // there can only be one segment that wraps, and that should be the last one (see sort_and_deoverlap() impl.)
                                        return true;
                                }
                        }

                        for (int32_t l{0}; l <= h;)
                        {
                                const auto m = (l & h) + ((l ^ h) >> 1);
                                const auto segment = segments[m];

                                if (segment.contains(token))
                                        return true;
                                else if (token <= segment.left)
                                        h = m - 1;
                                else
                                        l = m + 1;
                        }

                        return false;
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

                Ring(const std::vector<T> &v)
                    : Ring{v.data(), v.size()}
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

                const T &token_predecessor_by_index(const uint32_t idx) const noexcept
                {
                        return tokens[(idx + (cnt - 1)) % cnt];
                }

                const T &token_predecessor(const T token) const noexcept
                {
                        return token_predecessor_by_index(index_of(token));
                }

                const T &token_successor_by_index(const uint32_t idx) const noexcept
                {
                        return tokens[(idx + 1) % cnt];
                }

                const T &token_successor(const T token) const noexcept
                {
                        return token_successor_by_index(index_of(token));
                }

                auto index_segment(const uint32_t idx) const noexcept
                {
                        return ring_segment<T>(tokens[(idx + (cnt - 1)) % cnt], tokens[idx]);
                }

                // segment in the ring that owns this token
                // based on the (prev segment.token, this segment.token] ownership rule
                auto token_segment(const T token) const noexcept
                {
                        return index_segment(index_of(token));
                }

                // see also sort_and_deoverlap()
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

                auto tokens_segments(const std::vector<token_t> &t) const
                {
                        std::vector<segment_t> res;

                        res.reserve(t.size());
                        for (const auto token : t)
                        {
                                const auto idx = index_owner_of(token);

                                res.push_back({token_predecessor_by_index(idx), token});
                        }

                        std::sort(res.begin(), res.end(), [](const auto &a, const auto &b) { return a.left < b.left; });
                        return res;
                }

                // Assuming a node is a replica for tokens in segments `current`, and then it assumes ownership of a different
                // set of segments, `updated`
                //
                // This handy utility method will generate a pair of segments list:
                // 1. The first is segments the node will need to *fetch* from other nodes in the ring, because it will now be also responsible
                // for those segments, but it does not have the data, based on the current owned segments.
                // 2. The second is segments the node will need to *stream* to other nodes in the ring, because it will no longer hold data for them.
                //
                // Obviously, if a node is just introduced to a ring (i.e have only updated and no current segments ), it should
                // just fetch data for all the current segments. Conversely, if the node is exiting the ring, it should
                // consider streaming all the data it has to other nodes if needed, and not fetch any data to itself.
                //
                // make sure that current and updated are in order  e.g std::sort(start, end, [](const auto &a, const auto &b) { return a.left < b.left; });
                //
                // Because the output will be an array of segments (_not_ tokens), you will need to determine the segments of the ring that intersect it
                // in order to figure out which nodes have which parts of the segments.
                //
                // This is a fairly expensive method (although it should be easy to optimize it if necessary), but given how rare it should be used, that's not a real concern
                //
                // Example: current segment [10, 20), updated segment [10, 25)
                // Example: current segment [10, 20), updated segment [8, 30)
                static auto compute_segments_ownership_updates(const std::vector<segment_t> &currentSegmentsInput, const std::vector<segment_t> &updatedSegmentsInput)
                {
                        std::vector<segment_t> toFetch, toStream, current, updated, toFetchFinal, toStreamFinal;
                        segment_t segmentsList[2];

                        // We need to work on normalized lists of segments
                        current = currentSegmentsInput;
                        ring_segment<T>::sort_and_deoverlap(&current);

                        updated = updatedSegmentsInput;
                        ring_segment<T>::sort_and_deoverlap(&updated);

                        for (const auto curSegment : current)
                        {
                                const auto n = toStream.size();

                                for (const auto updatedSegment : updated)
                                {
                                        if (curSegment.intersects(updatedSegment))
                                                toStream.insert(toStream.end(), segmentsList, segmentsList + updatedSegment.difference(curSegment, segmentsList));
                                }

                                if (toStream.size() == n)
                                {
                                        // no intersection; accept whole segment
                                        toStream.push_back(curSegment);
                                }
                        }

                        for (const auto updatedSegment : updated)
                        {
                                const auto n = toFetch.size();

                                for (const auto curSegment : current)
                                {
                                        if (updatedSegment.intersects(curSegment))
                                                toFetch.insert(toFetch.end(), segmentsList, segmentsList + curSegment.difference(updatedSegment, segmentsList));
                                }

                                if (toFetch.size() == n)
                                {
                                        // no intersection; accept whole segment
                                        toFetch.push_back(updatedSegment);
                                }
                        }

                        // normalize output
                        ring_segment<T>::sort_and_deoverlap(&toFetch);
                        ring_segment<T>::sort_and_deoverlap(&toStream);

                        // mask segments:
                        // 	from segments to fetch, mask currently owned segments
                        //	from segments to stream, mask segments we will own (updated segments)
                        ring_segment<T>::mask_segments(toFetch, current, &toFetchFinal);
                        ring_segment<T>::mask_segments(toStream, updated, &toStreamFinal);

                        return std::make_pair(toFetchFinal, toStreamFinal);
                }

                // When a node acquires ring tokens(joins a cluster), it only disupts segments its token(s) fall into
                // Assuming a ring of tokens:  (10, 100, 150, 180, 200)
                // and a node joins a cluster, and acquires token 120
                // then it will only affect requests for (100, 120]
                // so it will need to fetch content for (100, 120] from somewhere. Where? well, from whichever owned (100, 150]
                // which is just the successor node, which we can find using index_owner_of()
                // This is a simple replication strategy implementation; we 'll just walk the ring clockwise and collect nodes that own
                // the tokens, skipping already collected nodes
                //
                // EXAMPLE: This is an illustrative example; you shouldn't really use this in production as is
                template <typename L>
                auto token_replicas_basic(const token_t token, const uint8_t replicationFactor, L &&endpoint_token) const
                {
                        using node_t = typename std::result_of<L(uint32_t)>::type;
                        std::vector<node_t> nodes;
                        const auto base = index_owner_of(token);
                        auto idx = base;

                        do
                        {
                                const auto node = endpoint_token(idx);

                                if (std::find(nodes.begin(), nodes.end(), node) == nodes.end())
                                {
                                        nodes.push_back(node);
                                        if (nodes.size() == replicationFactor)
                                                break;
                                }

                                idx = (idx + 1) % size();
                        } while (idx != base);

                        return nodes;
                }

                // This generates the lists of tokens and matching nodes that own them based on a new ownership state that results
                // from applying the changes in ringTokensNodes
                // Specifically, in the resulted topology, current nodes tokens are replaced with their updated set in ringTokensNodes
                template <typename node_t>
                std::pair<std::vector<token_t>, std::vector<node_t>> new_topology(const node_t *const ringTokensNodes,
                                                                                        const std::unordered_map<node_t, std::vector<token_t>> &futureNodesTokens) const
                {
                        std::vector<token_t> transientRingTokens;
                        std::vector<node_t> transientRingTokensNodes;
                        std::unordered_map<token_t, node_t> map;

                        for (uint32_t i{0}; i != cnt; ++i)
                        {
                                const auto token = tokens[i];

                                if (futureNodesTokens.find(ringTokensNodes[i]) == futureNodesTokens.end())
                                {
                                        transientRingTokens.push_back(tokens[i]);
                                        map.insert({tokens[i], ringTokensNodes[i]});
                                }
                        }

                        for (const auto &it : futureNodesTokens)
                        {
                                const auto node = it.first;

                                transientRingTokens.insert(transientRingTokens.end(), it.second.data(), it.second.data() + it.second.size());
                                for (const auto token : it.second)
                                        map.insert({token, node});
                        }

                        std::sort(transientRingTokens.begin(), transientRingTokens.end());

                        // The associated nodes for each token in the transient ring
                        transientRingTokensNodes.reserve(transientRingTokens.size());
                        for (const auto token : transientRingTokens)
                                transientRingTokensNodes.push_back(map[token]);

                        return {std::move(transientRingTokens), std::move(transientRingTokensNodes)};
                }

                template <typename node_t, typename L>
                static node_t *filter_by_distance(node_t *const nodes, const node_t *const end, L &&l)
                {
                        using dist_t = typename std::result_of<L(node_t)>::type;
                        dist_t min;
                        uint32_t out{0};

                        for (auto it = nodes; it != end; ++it)
                        {
                                if (!out)
                                {
                                        min = l(*it);
                                        nodes[out++] = *it;
                                }
                                else if (const auto d = l(*it); d == min)
                                        nodes[out++] = *it;
                                else if (d < min)
                                {
                                        min = d;
                                        nodes[0] = *it;
                                        out = 1;
                                }
                        }

                        return nodes + out;
                }


                // Whenever one more nodes alters the ring topology (when joining a cluster, leaving a cluster, or acquiring a different set of tokens), we need to
                // account for that change, by copying data to nodes that will serve segments they didn't already were serving(thus, they don't have the data for that ring space)
                // and by copying data to nodes that will now serve as a result of one or more nodes dropping segments they used to serve. This is necessary in order to
                // support replication semantics.
                //
                // You should initiate a transition, and when it is complete, create the final ring topology using
                // new_topology() like it's used here for the transient ring, and switch to it, by advertising the new tokens for all tokens in the ring.
                //
                // GUIDELINES
                // - There can be only one active transition in progress. If you allow for concurrent transitions, you will almost definitely end up with
                // invalid rings that likely contain missing data.
                // - For existing nodes participating in the transition: They should not be stopped or otherwise be treated in any special way.
                // - For nodes that are to join the cluster(i.e are not already participating in the ring), you should wait until the transition has completed successfully,
                // 	and then initialize them with the tokens you used for them in the transition.
                // - For nodes that are to leave the cluster(i.e notes in the current cluster, but not in the cluster topology after the transition), you should
                //	wait until the transition has completed successfully, and them stop them, and make sure you won't start them again with the same tokens.
                //
                // With this method, the only tricky operation becomes the coordindation required for switching to the new topology after the
                // streaming required for the transition is complete. You will need to (re)start nodes using specific tokens, etc.
                //
                // OPTIMIZATION OPPORTUNITIES
                // - You should use filter_by_distance() to filter sources, or a similar function so that you will always select among the closest(in terms of network hops) nodes to
                // the ring target node for streaming, in order to minimize streaming time.
                // - You should try to schedule the streaming operations fairly among the involved nodes. If you over-load a node and under-load the rest, or vice versa, the time
                // and effort(cost) will be much higher.
                //
                // EXAMPLES
                // - If you 'd like to add 5 new nodes to your cluster, you can pick appropriate tokens(functionality for selecting tokens from the ring based on current distribution will be
                // implemented later) for those new tokens, initiate a transition that involves them and their new tokens, and when you are done streaming, you should start those new 5 nodes, and each
                // should be conigured to use the tokens you selected for transition().
                // - If you 'd like to decomission 1 node, you just need to a new transition that involves that node, and the list of tokens it will own will be empty. Once the streaming is
                // complete, you should stop the node.
		//
		// MISC
		// Maybe it's a good idea to get all nodes accept all writes for the segments they serve, even if the transition hasn't been completed, so that you won't need
		// read-repair to deal with missing content. While the streaming is in-progress, writes for segments served by involved segments that they didn't already serve
		// in the current ring, will be lost - unless your service can support a mode where all writes eventually will be forwarded to those nodes, either by using some
		// sort of hinted-handoff system, or by having those nodes in the ring but only responding for reads for the current segments, and writes should get a diff. strategy.
		// This section will be eventually revised with more concrete advise on how to do about this.
                template <typename node_t, typename L>
                auto transition(
                    const node_t *const ringTokensNodes,
                    const std::unordered_map<node_t, std::vector<token_t>> &futureNodesTokens,
                    L &&replicas_for) const
                {
                        static constexpr size_t maxReplicasCnt{16};
                        const auto segments_of = [&replicas_for](const Ring &ring, const node_t *const ringTokensNodes, const node_t node, std::vector<segment_t> *const res) {
                                node_t replicas[maxReplicasCnt];

                                for (uint32_t i{0}; i != ring.cnt; ++i)
                                {
                                        const auto token = ring.tokens[i];
                                        const auto replicasCnt = replicas_for(ring, ringTokensNodes, token, replicas);

                                        if (std::find(replicas, replicas + replicasCnt, node) != replicas + replicasCnt)
                                                res->push_back({ring.token_predecessor_by_index(i), token});
                                }

                                std::sort(res->begin(), res->end(), [](const auto &a, const auto &b) { return a.left < b.left; });
                        };

                        const auto transientRingTopology = new_topology(ringTokensNodes, futureNodesTokens);
                        const auto &transientRingTokens = transientRingTopology.first;
                        const auto &transientRingTokensNodes = transientRingTopology.second;
                        const Ring transientRing(transientRingTokens.data(), transientRingTokens.size());
                        const auto transientRingSegments = transientRing.segments();
                        const auto currentRingSegments = segments();
                        std::vector<segment_t> outSegments;
                        segment_t segmentsList[2];
                        std::vector<std::pair<segment_t, std::pair<node_t, std::vector<node_t>>>> transportMap;
                        std::unordered_map<node_t, std::vector<segment_t>> curRingServeMap;
                        std::vector<node_t> replicas;
                        node_t tokenReplicas[maxReplicasCnt], futureReplicas[maxReplicasCnt];
			std::vector<segment_t> replicaForSegmentsFuture, replicaForSegmentsNow;

                        // Build (node => [segments]) map for the current ring
                        {
                                std::vector<std::pair<node_t, segment_t>> v;

                                for (const auto segment : currentRingSegments)
                                {
                                        const auto n = replicas_for(*this, ringTokensNodes, segment.right, tokenReplicas);

                                        for (uint8_t i{0}; i != n; ++i)
                                                v.push_back({tokenReplicas[i], segment});
                                }

                                std::sort(v.begin(), v.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

                                const auto n = v.size();
                                const auto all = v.data();

                                for (uint32_t i{0}; i != n;)
                                {
                                        const auto node = v[i].first;
                                        std::vector<segment_t> list;

                                        do
                                        {
                                                list.push_back(v[i].second);
                                        } while (++i != n && v[i].first == node);

                                        curRingServeMap.insert({node, std::move(list)});
                                }
                        }

                        for (const auto &it : futureNodesTokens)
                        {
                                const auto node = it.first;

				replicaForSegmentsFuture.clear();
				replicaForSegmentsNow.clear();

                                segments_of(transientRing, transientRingTokensNodes.data(), node, &replicaForSegmentsFuture);
                                segments_of(*this, ringTokensNodes, node, &replicaForSegmentsNow);

                                // Compute what needs to be delivered to _this_ node
                                for (const auto futureSegment : replicaForSegmentsFuture)
                                {
                                        //  Mask segments this node serves already, no need to acquire any content we already have
					const auto futureSegmentWraps = futureSegment.wraps();
                                        outSegments.clear();
                                        segment_t::mask_segments(&futureSegment, (&futureSegment) + 1, replicaForSegmentsNow, &outSegments);


                                        if (outSegments.empty())
                                        {
                                                // No need to acquire extra data
                                                continue;
                                        }

					// TODO: use binary search to locate the next segment in currentRingSegments
					// No need for the linear scan


                                        // OK, so who's going to provide content for those segments, based on the current ring?
                                        for (const auto it : currentRingSegments)
                                        {
						if (it.right <= futureSegment.left)
						{
							// can safely skip it
							continue;
						}
						else if (it.left > futureSegment.right && !futureSegmentWraps && !it.wraps())
						{
							// can safely stop here
							break;
						}

                                                const auto cnt = it.intersection(futureSegment, segmentsList);

                                                if (!cnt)
                                                        continue;

                                                const std::vector<node_t> replicas(tokenReplicas, tokenReplicas + replicas_for(*this, ringTokensNodes, it.right, tokenReplicas));

                                                for (uint8_t i{0}; i != cnt; ++i)
                                                        transportMap.push_back({segmentsList[i], {node, replicas}}); // from replicas, to node, that segment
                                        }
                                }

                                // Whenever a node gives up (part of) a ring segment, we need to shift data around in order
                                // to account for the fact that replication factor for the data that span that segment will drop by 1.
                                for (const auto currentSegment : replicaForSegmentsNow)
                                {
                                        const auto token = currentSegment.right;
					const auto currentSegmentWraps = currentSegment.wraps();
                                        bool haveSources{false};

					// TODO: use binary search to locate the next segment in transientRingSegments
					// No need for linear scan

                                        for (const auto futureSegment : transientRingSegments)
                                        {
						if (futureSegment.right <= currentSegment.left)
						{
							// can safely skip it
							continue;
						}
						else if (futureSegment.left > currentSegment.right && !currentSegmentWraps && !futureSegment.wraps())
						{
							// can safely stop here
							break;
						}

                                                const auto cnt = futureSegment.intersection(currentSegment, segmentsList);

                                                if (!cnt)
                                                        continue;

                                                const auto futureReplicasCnt = std::remove_if(futureReplicas,
                                                                                              futureReplicas + replicas_for(transientRing, transientRingTokensNodes.data(), futureSegment.right, futureReplicas),
                                                                                              [node, &futureNodesTokens](const node_t target) {
                                                                                                      if (target == node)
                                                                                                      {
                                                                                                              // exclude self
                                                                                                              return true;
                                                                                                      }
                                                                                                      else if (futureNodesTokens.find(target) != futureNodesTokens.end())
                                                                                                      {
                                                                                                              // exclude nodes that are also involved in this process, otherwise we may output the same value twice in transportMap
                                                                                                              return true;
                                                                                                      }
                                                                                                      else
                                                                                                      {
                                                                                                              return false;
                                                                                                      }

                                                                                              }) -
                                                                               futureReplicas;

                                                for (uint8_t i{0}; i != cnt; ++i)
                                                {
                                                        const auto subSegment = segmentsList[i]; // intersection

							for (uint32_t ri{0}; ri != futureReplicasCnt; ++ri)
                                                        {
                                                                const auto target = futureReplicas[ri];

                                                                if (!haveSources)
                                                                {
                                                                        // lazy generation of the sources(replicas) for this segment
                                                                        // replicas should include this node
                                                                        replicas.clear();
                                                                        replicas.insert(replicas.end(), tokenReplicas, tokenReplicas + replicas_for(*this, ringTokensNodes, token, tokenReplicas));
                                                                        haveSources = true;
                                                                }

                                                                if (auto s = curRingServeMap.find(target); s != curRingServeMap.end())
                                                                {
                                                                        // this node serves 1+ segments already in the current segment
                                                                        // mask subSegment with them; we don't want to send data to nodes if they already have any of it
                                                                        outSegments.clear();
                                                                        segment_t::mask_segments(&subSegment, (&subSegment) + 1, s->second, &outSegments);

                                                                        for (const auto s : outSegments)
                                                                                transportMap.push_back({s, {target, replicas}});
                                                                }
                                                                else
                                                                {
                                                                        // this target does not currently server any segments in the current segment
                                                                        transportMap.push_back({subSegment, {target, replicas}});
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }

                        return transportMap;
                }
        };
}

#ifdef HAVE_SWITCH
template <typename token_t>
static inline void PrintImpl(Buffer &b, const ConsistentHashing::ring_segment<token_t> &segment)
{
        b.append("(", segment.left, ", ", segment.right, "]");
}

template <typename T>
static inline void PrintImpl(Buffer &b, const ConsistentHashing::Ring<T> &ring)
{
        b.append(_S32("(( "));
        if (const auto cnt = ring.cnt)
        {
                for (uint32_t i{1}; i != cnt; ++i)
                        b.append(ConsistentHashing::ring_segment<T>(ring.tokens[i - 1], ring.tokens[i]), ",");

                b.append(ConsistentHashing::ring_segment<T>(ring.tokens[cnt - 1], ring.tokens[0]));
        }
        b.append(_S32(" ))"));
}
#endif
