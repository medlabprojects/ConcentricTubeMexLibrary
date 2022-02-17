/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#ifdef WIN32
#include "Windows.h"
#endif

#ifdef __GLIBC__
#include <time.h>
#endif

namespace Utility
{

	class CodeTimer
	{
	public:
		CodeTimer()
		{
#ifdef WIN32
			QueryPerformanceFrequency(&freq);
#endif
#ifdef __GLIBC__
			clock_getres(CLOCK_REALTIME, &freq);	
#endif
		}

		void start() {
#ifdef WIN32
			QueryPerformanceCounter(&t1);
#endif
#ifdef __GLIBC__
			clock_gettime(CLOCK_REALTIME, &t1);
#endif
		}

		void stop() {
#ifdef WIN32
			QueryPerformanceCounter(&t2);
#endif
#ifdef __GLIBC__
			clock_gettime(CLOCK_REALTIME, &t2);
#endif
		}

		double elapsed() {
#ifdef WIN32
			return static_cast<double>(t2.QuadPart - t1.QuadPart) / static_cast<double>(freq.QuadPart);
#endif
#ifdef __GLIBC__
			return static_cast<double>(t2.tv_sec - t1.tv_sec) + 
						static_cast<double>(t2.tv_nsec - t1.tv_nsec)/1e9;
#endif
		}

	private:
#ifdef WIN32
		LARGE_INTEGER t1, t2, freq;
#endif
#ifdef __GLIBC__
		timespec freq, t1, t2;
#endif
	};


}
