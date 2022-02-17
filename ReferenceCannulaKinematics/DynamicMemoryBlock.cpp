#include "DynamicMemoryBlock.h"

GlobalMemoryBlock *pGlobalMemoryBlock__;

void InitGlobalMemoryBlock()
{
	pGlobalMemoryBlock__ = new GlobalMemoryBlock();
}

void FinalizeGlobalMemoryBlock()
{
	delete pGlobalMemoryBlock__;
}

void ResetGlobalMemoryBlock()
{
	pGlobalMemoryBlock__->reset();
}