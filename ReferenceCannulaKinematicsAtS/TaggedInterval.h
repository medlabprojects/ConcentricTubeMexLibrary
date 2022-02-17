/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <tuple>
#include <vector>
#include <algorithm>
#include <boost/container/static_vector.hpp>

/**
	\brief Tagged Interval handling for generating integration points along the robot length

	This namespace contains all of the functions needed to analyze the integration waypoints
	for a concentric tube robot. The primary routine is TInterval::resolve_to_disjoint, which takes a list
	of intervals of the form [a,b,tag], where [a,b] is an interval on the real line and the 
	tag indicates a desired point density inside this interval, and returns a list of
	"disjoint" or abutting intervals so that each point is assigned the highest density
	needed. As the tubes of a concentric tube robot translate with respect to one another, the
	natural integration waypoints defined by each tube translate in the arc length space, 
	and this namespace provides the necessary tools to correctly compute integration waypoints.

	An example and a corresponding code listing are shown in the documentation of TInterval::resolve_to_disjoint.
*/
namespace TInterval {
	const static int MAX_INTERVALS = 50;

	/**
		\brief Enumerated point densities
	*/
	typedef enum {
		SPARSE = 0,
		NORMAL = 1,
		DENSE = 2,
		DISCONTINUITY = 3
	} IntervalDensity;

	/**
		\brief Interval endpoint names
	*/
	typedef enum {
		LEFT = 0,
		RIGHT
	} IntervalEnd;

	/**
		\brief An Interval on the real line together with a tag

		This typedef is an alias for a tuple that holds the two endpoint intervals
		and the tag or IntervalDensity. The three components can be retrieved with
		either std::get or with the convenience functions left(TaggedInterval),
		right(TaggedInterval), and density(TaggedInterval).
	*/
	typedef std::tuple<double, double, IntervalDensity>	TaggedInterval; 

	/**
		\brief An interval endpoint with a tag

		The ends of intervals are identified by three items: the location,
		whether it is the left or right end of the interval, and the tag associated
		with the interval.
	*/
	typedef std::tuple<double, IntervalEnd, IntervalDensity> TaggedEndpoint; 

	/**
		\brief A memory block

		This memory block is a statically allocated array of 8192 bytes,
		which is used by the interval arithmetic code as scratch space
		for dynamically sized vectors (no calls to malloc/free take place in namespace TInterval).
	*/
	//typedef StackMemoryBlock<8192>	MemBlock;

	/**
		\brief A Stack allocator for the std containers, such as vector

		This allocator provides the standard stl containers with memory
		out of a MemBlock. Example: 

		~~~{.cpp}
			TInterval::MemBlock mb;
			TInterval::IntervalAlloc alloc(mb);
			TInterval::IntervalList interval_list(alloc);
		~~~

	*/
	//typedef ::StackAllocator< TaggedInterval >	IntervalAlloc;
	//typedef ::StackAllocator< TaggedEndpoint >    EndpointAlloc;
	//typedef ::StackAllocator< IntervalDensity >   TagAlloc;
	/**
		\brief A list of tagged intervals

		Lists of tagged intervals are held in a standard vector container, which
		has been modified to use a stack memory block thorugh the ::StackAllocator class.
	*/
	//typedef std::vector < TaggedInterval, IntervalAlloc > IntervalList;
	typedef boost::container::static_vector< TaggedInterval, MAX_INTERVALS >	IntervalList;

	/**
		\brief A list of tagged endpoints

		Lists of tagged endpoints are held in a standard vector container, which
		has been modified to use a stack memory block thorugh the ::StackAllocator class.
	*/
	//typedef std::vector < TaggedEndpoint, EndpointAlloc >	EndpointList;
	typedef boost::container::static_vector< TaggedEndpoint, 2*MAX_INTERVALS >	EndpointList;

	/**
		\brief A heap structure of tagged endpoints

		This heap is really just a standard vector, but is typedef'ed for clarity.
	*/
	//typedef std::vector < IntervalDensity, TagAlloc > TagHeap;
	typedef boost::container::static_vector< IntervalDensity, 2*MAX_INTERVALS >	TagHeap;

	/**
		\brief Get the left endpoint location of a tagged interval
		\param t a tagged interval
		\return the left endpoint of the interval
	*/
	double left(TaggedInterval const& t);

	/**
		\brief Get the right endpoint location of a tagged interval
		\param t a tagged interval
		\return the right endpoint of the interval
	*/
	double right(TaggedInterval const& t);

	/**
		\brief Get the length of a tagged interval
		\param t a tagged interval
		\return the length of the interval

		This function is a convenience method for right(t)-left(t)
	*/
	double length(TaggedInterval const& t);

	/**
		\brief Get the density or "tag" of a tagged interval
		\param t a tagged interval
		\return the density tag of the interval
	*/
	IntervalDensity density(TaggedInterval const& t);

	/**
		\brief Get the entire left endpoint of the interval as a TaggedEndpoint
		\param t a tagged interval
		\return an endpoint at the left end of the interval with the same tag
	*/
	TaggedEndpoint left_end(TaggedInterval const& t);

	/**
		\brief Get the entire right endpoint of the interval as a TaggedEndpoint
		\param t a tagged interval
		\return an endpoint at the right end of the interval with the same tag
	*/
	TaggedEndpoint right_end(TaggedInterval const& t);

	/**
		\brief Get the location of an endpoint
		\param e a tagged endpoint
		\return the location of a tagged endpoint
	*/
	double endpoint_location(TaggedEndpoint const& e);

	/**
		\brief Get the type (LEFT, RIGHT) of an endpoint
		\param e a tagged endpoint
		\return whether the endpoint is a LEFT endpoint or a RIGHT endpoint
	*/
	IntervalEnd endpoint_type(TaggedEndpoint const& e);

	/**
		\brief Get the "tag" or density of an endpoint
		\param e a tagged endpoint
		\return the density associated with a tagged endpoint
	*/
	IntervalDensity endpoint_density(TaggedEndpoint const& e);

	/**
		\brief Compare two endpoints
		\param t1 The first endpoint
		\param t2 The second endpoint
		This function compares first the location, and then the
		type. The tag does not participate in the comparison.
	*/
	bool operator<(TaggedEndpoint const& t1, TaggedEndpoint const& t2);

	/**
		\brief Convert a list of intervals to a list of endpoints
		\param IL the interval list
		This "splits" the endpoints of intervals apart and packs them
		all into an EndpointList. This function is used in the process
		of converting intervals to a disjoint list of intervals with
		resolve_to_disjoint
	*/
	EndpointList convert_to_endpoints(IntervalList const& IL); 

	/**
		\brief Get the total length of a list of intervals (assumed to be disjoint intervals)
		\param IL	the interval list (disjoint)
		\return the total length of the list of intervals
	*/
	double total_length(IntervalList const& IL);

	/**
		\brief	Adds a tag to the TagHeap (called for Left endpoints)
	*/
	bool add_tag(TagHeap& Tags, IntervalDensity tag, IntervalDensity &old_max);

	/**
		\brief	Removes a tag to the TagHeap (called for Right endpoints)
	*/
	bool remove_tag(TagHeap& Tags, IntervalDensity tag, IntervalDensity &old_max);

	/**
		\brief	Accepts endpoints and dispatches to add_tag and remove_tag
	*/
	bool AcceptEndpoint(TagHeap& Tags, TaggedEndpoint const& E, IntervalDensity& old_max);

	/**
		\brief	Resolves a list of TaggedIntervals to a disjoint list
		\param IL The IntervalList to resolve
		\return The disjoint list of intervals, prioritized by tag

		This function is the main routine of the TaggedInterval component of the
		library. Resolving the IntervalList to disjoint intervals results in a list
		of intervals with abutting endpoints, where each element of the real line
		in between the left-most endpoint in IL and the right-most endpoint in IL is
		assigned the highest tag (most density) that is associated with that point
		in any of the input intervals. As an example, consider the intervals

		[0,1,DENSE], [0.5,2,NORMAL], [1.5,3,DENSE], [2,10,SPARSE]

		These intervals must be resolved so that for each point between 0 and 10, 
		the point belongs to an interval with the highest density that it has in
		any of the input intervals. The correct answer is given by the list

		[0,1,DENSE], [1,1.5,NORMAL], [1.5,3,DENSE], [3,10,SPARSE]

		Complete example code corresponding to this example, which
		will print the result to cout:

		~~~{.cpp}
			TInterval::TaggedInterval ti1 (0.0, 1.0, TInterval::DENSE);
			TInterval::TaggedInterval ti2 (0.5, 2, TInterval::NORMAL);
			TInterval::TaggedInterval ti3 (1.5, 3.0, TInterval::DENSE);
			TInterval::TaggedInterval ti4 (2.0, 10.0, TInterval::SPARSE);
			TInterval::MemBlock mb;
			TInterval::IntervalAlloc alloc(mb);
			TInterval::IntervalList IL(alloc);
			IL.push_back(ti1); 
			IL.push_back(ti2); 
			IL.push_back(ti3);
			IL.push_back(ti4);
			TInterval::IntervalList DL = TInterval::resolve_to_disjoint(IL);
			std::cout << DL << std::endl;
		~~~
	*/
	IntervalList resolve_to_disjoint(IntervalList const& IL); 

	/**
		\brief Stream insertion for TaggedInterval objects (text)
	*/
	template <typename Stream>
	Stream& operator<<(Stream& os, TInterval::TaggedInterval const& I)
	{
		os << "(" << TInterval::left(I) << "," << TInterval::right(I) << "," << TInterval::density(I) << ")";
		return os;
	}

	/**
		\brief Stream insertion for an IntervalList object (text)
	*/
	template <typename Stream>
	Stream& operator<<(Stream& os, TInterval::IntervalList const& IL)
	{
		std::copy(std::begin(IL), std::end(IL), std::ostream_iterator<TInterval::TaggedInterval>(os, " : \n"));
		return os;
	}

	struct has_tag
	{
		has_tag(IntervalDensity d_) : d(d_) {}

		bool operator()(TaggedInterval const& ival)
		{
			if (density(ival) == d) return true;
			else return false;
		}
	private:
		IntervalDensity	d;
	};

	IntervalList	intersect(IntervalList const& IL, TInterval::TaggedInterval const& ival);
}


