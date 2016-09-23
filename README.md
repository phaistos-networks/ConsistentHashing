This is an C++14 implementation of [Consistent hashing](https://en.wikipedia.org/wiki/Consistent_hashing), abstracted as a `Ring` of tokens, with a `ring_segment` data structure that represents a segment of the that ring. [We](http://phaistosnetworks.gr/) have been using this implementation for many years on multiple distributed systems, including our massive scale, high performance distributed store (CloudDS). 

Please check the comments in the single header file for how to use it and how it works. It is pretty trivial to use it and various useful methods are implemented for building robust distributed services.

