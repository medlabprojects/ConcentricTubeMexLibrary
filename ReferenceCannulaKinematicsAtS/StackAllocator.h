/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once
#include <cstdint>
#include <cstddef>
#include <iterator>
#include <memory>
#include "StackMemoryBlock.h"
#include "DynamicMemoryBlock.h"
#include <functional>

/**
	\brief The StackAllocator class is a stack-memory based standard allocator that uses the global memory provided by GlobalMemoryBlock

	This class allows the standard containers to be used without fear of calls
	to malloc/free, which are typically not allowed in hard real-time code.

	\param T The type to be allocated, usually matches the container type but does not have to.
	\param MemBlockT The memory block, which should usually be StackMemoryBlock<N> for N bytes
			of storage.

	Note that this allocator is designed for each area of memory to be "consumed" as it is 
	used, so vector expansions, which take place as new allocations, can rapidly use up
	all of the available memory. This allows the allocator to be simple and fast, but 
	is a tradeoff requiring increased memory storage space. On modern processors with
	many megabytes of cache space, this is unlikely to be an issue for most of the tasks
	performed by this allocator. If the allocator is used with the standard vector class
	to store a large amount of data, it would be wise to make an appropriate call to .reserve() 
	for the container using the	allocator, to avoid automatic resizing.
	*/
template <typename T>
class StackAllocator
{
public:
	typedef T			value_type;
	typedef T*			pointer;
	typedef T&			reference;
	typedef const T*	const_pointer;
	typedef const T&	const_reference;
	typedef size_t		size_type;
	typedef ptrdiff_t	difference_type;

	template <typename W> struct rebind {
		typedef StackAllocator<W> other;
	};

	//! Required by STL
	pointer address(reference x) const {
		return &x;
	}

	//! Required by STL
	const_pointer address(const_reference x) const {
		return &x;
	}

	//! Allocation is simple, and just gives the next chunk of available memory
	//The hint is always ignored
	pointer allocate(size_type n, std::allocator<void>::const_pointer hint = 0)
	{
		return reinterpret_cast<T*>(pGlobalMemoryBlock__->template allocate< value_type >(n));
	}

	//! No deallocation ever takes place in the StackAllocator, but this is required by STL
	void deallocate(pointer p, size_type n)
	{}

	//! Returns the maximum number of storage bytes
	size_type max_size() const
	{
		return pGlobalMemoryBlock__->size();
	}

	//! Calls placement new to copy construct a new element in the storage space
	template <class U>
	void construct(U* p, U const& val)
	{
		new ((void*)p) U(val);
	}

	//! Calls the object destructor
	template <class U>
	void destroy(U* p)
	{
		p->~U();
	}

	template <typename U>
	friend class StackAllocator ;
	
	/**
		\brief Default constructor
	*/
	StackAllocator() {}

	/**
		\brief Copy construct a StackAllocator from another StackAllocator, using the same storage space

		This allows a container for any type to use an allocator for any other type. (STL requirement)
	*/
	template <typename U>
	StackAllocator(const StackAllocator<U>& alloc)
	{}

private:

	template <typename T1, typename T2>
	friend bool operator==(StackAllocator<T1> const& alloc1, StackAllocator<T2> const& alloc2);
};

template <typename T1, typename T2>
bool operator==(StackAllocator<T1> const& alloc1, StackAllocator<T2> const& alloc2)
{
	return true;
}


template <typename T1, typename T2>
bool operator!=(StackAllocator<T1> const& alloc1, StackAllocator<T2> const& alloc2)
{
	return !(alloc1 == alloc2);
}