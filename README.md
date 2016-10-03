This is an C++14 implementation of [Consistent hashing](https://en.wikipedia.org/wiki/Consistent_hashing), abstracted as a `Ring` of tokens, with a `ring_segment` data structure that represents a segment of the ring. [We](http://phaistosnetworks.gr/) have been using this implementation for many years in building multiple distributed systems, including our massive scale, high performance distributed store (CloudDS). 

Please check the comments in the single header file for how to use the data structures and their APIs, and how it works. 
It is pretty trivial to use it and various useful methods are implemented for building robust distributed services.

### Using it in your Project
Just include [consistent_hashing.h](https://github.com/phaistos-networks/ConsistentHashing/blob/master/consistent_hashing.h), and make sure you set `std=c++14` or higher compiler option.

If you are going to need more than 2^64 ring tokens, and you probably should, you will need a struct or class to represent it (because an uin64_t won't suffice). In this case, you need to implement a few things:

- You need a `TrivialCmp()` implementation for your token type. This should return < 0, 0, or > 0 depending on the comparison result of two tokens. Hopefully, a future C++ standard update will introduce the [spaceship operator](https://en.wikipedia.org/wiki/Three-way_comparison) we could override and solve this more elegantly, but for now this will have to do.
```cpp
template<>
static inline int8_t TrivialCmp<hugetoken_t>(const hugetoken_t &a, const hugetoken_t &b)
{
	// return comparison result
}
```

- You will need to implement an appropriate `std::numeric_limits<hugetoken_t>::min()`, like so:
```cpp
namespace std
{
	template<>
	struct numeric_limits<hugetoken_t>
	{
		static inline const hugetoken_t min()
		{
			//return minimum possible token (e.g 0)
		}
	};
}
```

- Implement appropriate `std::min()` and `std::max()` methods like so:
```cpp
namespace std
{
	template<>
	inline const hugetoken_t &min<hugetoken_t>(const hugetoken_t &a, const hugetoken_t &b)
	{
		// return whichever is lower
	}

	template<>
	inline const hugetoken_t &max<hugetoken_t>(const hugetoken_t &a, const hugetoken_t &b)
	{
		// return whichever is higher 
	}
}
```

- Finally, you should implement appropriate `operator==`, `operator!=`, `operator<` and `operator>` methods for your hugetoken_t

This is really not that much work, and chances are you are already doing that anyway to support other needs of your codebase.

### Example

```cpp
#include <consistent_hashing.h>

int main()
{
        using token_t = uint32_t;
        using segment_t = ConsistentHashing::ring_segment<token_t>;
        using ring_t = ConsistentHashing::Ring<token_t>;
        using node_t = uint32_t;
        // Suppose we have a simple ring, and of those tokens, only
        // one is owned by the node we wish to update (node 1)
        std::pair<node_t, token_t> ringStructure[] =
            {
                {100, 10},
                {200, 20},
                {300, 30},
                {400, 40},
                {500, 50},
                {600, 60},
                {1, 70}, /* this is the only token owned by node 1 */
                {800, 80},
                {900, 90},
                {1000, 100},
                {150, 110},
                {1200, 120}};
        std::vector<node_t> ringTokensNodes, ringTokens;

        // We need all ring tokens, and all associated nodes for those tokens.
        // We also need to collect the tokens owned by the node
        for (auto it = std::begin(ringStructure), end = std::end(ringStructure); it != end; ++it)
        {
                const auto &v = *it;

                ringTokens.push_back(v.second);
                ringTokensNodes.push_back(v.first);
        }

        // This is a simple method that returns the replicas for a given token
        // In practice, you wouldn't allocate any memory here (i.e no std::vector<> use), you 'd
        // pay a lot of attention to performance, and you 'd possiblyh consider physical placement of node
        // (i.e at least one from local DC and the rest from other DCs)
        const auto replicas_of = [](const auto &ring, const auto ringTokensNodes, const token_t token, node_t *const out) {
                static constexpr uint8_t replicationFactor{2}; // how many copies of each ring segment we need at any given time
                const auto base = ring.index_owner_of(token);  // index in the ring tokens for the token-owner of token(input)
                uint32_t i{base}, n{0};

                // walk the the ring clockwise until we have enough nodes to return
                do
                {
                        const auto node = ringTokensNodes[i];

                        if (std::find(out, out + n, node) == out + n)
                        {
                                // haven't collected that node yet
                                // we only care for distinct nodes
                                out[n++] = node;
                                if (n == replicationFactor)
                                        break;
                        }

                        i = (i + 1) % ring.size();
                } while (i != base);

                return n;
        };

        // This is our ring
        const ring_t ring(ringTokens.data(), ringTokens.size());
        // those are the changes
        const std::unordered_map<node_t, std::vector<token_t>> topologyUpdates{
            {1, {35, 95}},
	    {110, {}}, 	// will remove node from the ring
            {64, {7}}};
        // Figure out what needs to be transfered, the available sources for those segments, and the targets
        auto updates = ring.transition(ringTokensNodes.data(), topologyUpdates, replicas_of);


	for (auto &it : updates)
        {
                const auto segment = it.first;
                const auto &toFrom = it.second;
                const auto target = toFrom.first;
                auto sources = toFrom.second;

                // Let's filter the sources, so that we only pull from the nodes closest to us in terms of node hopes
                sources.resize(ring_t::filter_by_distance(sources.data(), sources.data() + sources.size(), [target](const auto node) {
                                       return 1; // TODO: return an appropriate distance from target to node
                               }) -
                               sources.data());
        }

        // This function should consider all sources, consider opportunities for fairly scheduling of transfers
        // across the distinct sources, and otherwise do what it takes to transfer data among nodes
        //
        // One way to do this is to register a new "transfer session", and coordinate the process among all
        // ring nodes involved. In case of failure, the participating nodes should also abort.
        // You should also make sure that only one transfer is active at any given time.
        //
        // Eventually, the callback should be invoked, and that callback
        // should create a new final ring

        schedule_transfer(std::move(updates),
                          [ ringTokensNodes, ring, topologyUpdates ]() {
                                  const auto newTopology = ring.new_topology(ringTokensNodes.data(), topologyUpdates);

                                  switch_to_ring(newTopology);
                          });
        return 0;
}
```
--

With the included data structures and implemented algorithms, it should be trivial to build robust consistent-hashing based replication for your distributed systems.

Have Fun!
