/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <cstdint>
#include <cstdlib>
#include <stdexcept>

#ifdef WIN32
#define ALIGN __declspec(align(16))
#else
#define ALIGN
#endif

/**
\brief This class acts as a stack-allocated memory resourse for the StackAllocator allocator class

\param N The number of kilobytes of storage to allocate

Note: If an allocation is attempted which would result in a memory access violation,
std::bad_alloc is thrown.
*/
template <size_t N, size_t Align = 1>
class DynamicMemoryBlock
{
public:
	typedef uint8_t value_type;
	typedef uint8_t* pointer;
	typedef uint8_t& reference;
	typedef const pointer const_pointer;
	typedef const reference const_reference;
	typedef size_t	size_type;
	typedef ptrdiff_t	difference_type;

	DynamicMemoryBlock() {
		m_Mem = (pointer)malloc(1024*N + Align);
		void *aligned_base = reinterpret_cast<void*>(reinterpret_cast<size_t>(m_Mem) & ~(size_t(Align-1)));
		m_P = reinterpret_cast<pointer>(aligned_base) + (Align == 1 ? 0 : Align);

		if (m_Mem == nullptr)
		{
			throw std::bad_alloc();
		}
	}
	~DynamicMemoryBlock() {
		free(m_Mem);
	}

	pointer& GetP() { return m_P; }
	const_pointer GetM() { return m_Mem; }
	size_type size() { return 1024*N; }
	const_pointer end() { return m_Mem + size(); }

	inline bool check_alloc(size_type n) {
		if (m_P + n > m_Mem + size())
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

		if (bytes % Align != 0) { //need to realign the pointer
			size_type offset = Align - (bytes % Align);
			m_P += offset;
		}
		return p;
	}

	void reset()
	{
		void *aligned_base = reinterpret_cast<void*>(reinterpret_cast<size_t>(m_Mem)& ~(size_t(Align - 1)));
		m_P = reinterpret_cast<pointer>(aligned_base)+(Align == 1 ? 0 : Align);
	}
private:
	pointer m_Mem;
	pointer m_P;
};

typedef DynamicMemoryBlock<4096, 16>	GlobalMemoryBlock;
extern GlobalMemoryBlock	*pGlobalMemoryBlock__;

void InitGlobalMemoryBlock();
void FinalizeGlobalMemoryBlock();
void ResetGlobalMemoryBlock();