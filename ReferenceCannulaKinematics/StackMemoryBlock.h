/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <cstdint>

#ifdef WIN32
#define ALIGN __declspec(align(16))
#else
#define ALIGN
#endif

/**
	\brief This class acts as a stack-allocated memory resourse for the StackAllocator allocator class

	\param N The number of bytes of storage to allocate

	Note: If an allocation is attempted which would result in a memory access violation,
	std::bad_alloc is thrown.
*/
template <size_t N>
class StackMemoryBlock
{
public:
	typedef unsigned char value_type;
	typedef unsigned char* pointer;
	typedef unsigned char& reference;
	typedef const unsigned char* const_pointer;
	typedef const unsigned char& const_reference;
	typedef size_t	size_type;
	typedef ptrdiff_t	difference_type;

	StackMemoryBlock() : m_P(&(m_Mem[0])) {}
	~StackMemoryBlock() {}

	pointer& GetP() { return m_P;  }
	const_pointer GetM() { return m_Mem; }
	size_type size() { return N; }
	const_pointer end() { return std::end(m_Mem); }

	inline bool check_alloc(size_type n) {
		if (m_P + n >= std::end(m_Mem))
			return false;
		else
			return true;
	}
	
	template <typename T>
	pointer allocate(size_type n) {
		size_type bytes = n*sizeof(T);
		if (!check_alloc(bytes)) throw std::bad_alloc();
		pointer p = m_P;
		m_P += bytes;
		return p;
	}

private:
	ALIGN unsigned char m_Mem[N];
	pointer m_P;
};

