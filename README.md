This is an C++14 implementation of [Consistent hashing](https://en.wikipedia.org/wiki/Consistent_hashing), abstracted as a `Ring` of tokens, with a `ring_segment` data structure that represents a segment of the that ring. [We](http://phaistosnetworks.gr/) have been using this implementation for many years on multiple distributed systems, including our massive scale, high performance distributed store (CloudDS). 

Please check the comments in the single header file for how to use it and how it works. It is pretty trivial to use it and various useful methods are implemented for building robust distributed services.

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
};
```

- Finally, you should implement appropriate `operator==`, `operator!=`, `operator<` and `operator>` methods for your hugetoken_t

This is really not that much work, and chances are you are already doing that anyway to support other needs of your codebase.

--

With the included data structures and implemented algorithms, it should be trivial to build robust consistent-hashing based replication for your distributed systems.

Have Fun!
