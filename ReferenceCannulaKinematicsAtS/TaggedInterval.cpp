#include "TaggedInterval.h"
#include <numeric>
#include "DynamicMemoryBlock.h"
#include <functional>
#include <iterator>
#include <algorithm>
#include <iostream>

namespace TInterval {

	using std::get;
	using std::push_heap;
	using std::make_heap;

	double left(TaggedInterval const& t)
	{
		return get<0>(t);
	}

	double right(TaggedInterval const& t)
	{
		return get<1>(t);
	}

	IntervalDensity density(TaggedInterval const& t)
	{
		return get<2>(t);
	}

	TaggedEndpoint left_end(TaggedInterval const& t)
	{
		return TaggedEndpoint(left(t), LEFT, density(t));
	}

	TaggedEndpoint right_end(TaggedInterval const& t)
	{
		return TaggedEndpoint(right(t), RIGHT, density(t));
	}

	double endpoint_location(TaggedEndpoint const& e)
	{
		return get<0>(e);
	}

	IntervalEnd endpoint_type(TaggedEndpoint const& e)
	{
		return get<1>(e);
	}

	IntervalDensity endpoint_density(TaggedEndpoint const& e)
	{
		return get<2>(e);
	}

	bool operator<(TaggedEndpoint const& t1, TaggedEndpoint const& t2)
	{
		if (endpoint_location(t1) < endpoint_location(t2))
			return true;
		else if (endpoint_location(t1) == endpoint_location(t2) &&
			endpoint_type(t1) == LEFT && endpoint_type(t2) == RIGHT)
			return true;
		else if (endpoint_location(t1) == endpoint_location(t2) &&
			endpoint_type(t1) == endpoint_type(t2) &&
			endpoint_density(t1) < endpoint_density(t2))
			return true;
		else
			return false;
	}

	EndpointList convert_to_endpoints(IntervalList const& IL) {
		EndpointList EL;
		std::for_each(std::begin(IL), std::end(IL),
			[&EL](TaggedInterval const& I) {
			EL.push_back(left_end(I));
			EL.push_back(right_end(I));
		});

		std::sort(std::begin(EL), std::end(EL));
		return EL;
	}

	bool add_tag(TagHeap& Tags, IntervalDensity tag, IntervalDensity &old_max)
	{
		if (Tags.size() == 0) { //the queue is empty
			Tags.push_back(tag);
			return false;
		}
		else { //the queue is not empty
			old_max = Tags.front();
			Tags.push_back(tag);
			push_heap(Tags.begin(), Tags.end());
			if ( (tag > old_max) || 
              (tag == DENSE && old_max == DENSE)) //opening tag is higher or dense
				return true;
			else
				return false;
		}
	}

	bool remove_tag(TagHeap& Tags, IntervalDensity tag, IntervalDensity &old_max)
	{
		old_max = Tags.front();
		auto it = std::find(Tags.begin(), Tags.end(), tag);
		if (it != Tags.end())
			Tags.erase(it);

		if (Tags.size() > 0) {
			make_heap(Tags.begin(), Tags.end());
			unsigned char new_max = Tags.front();
			if (new_max < old_max || tag == DENSE) //we closed the highest tag (or it was dense)
				return true;
			else
				return false;
		}
		else { //we removed the last tag in the queue, so this closes an interval
			return true;
		}
	}

	bool AcceptEndpoint(TagHeap& Tags, TaggedEndpoint const& E, IntervalDensity& old_max)
	{
		IntervalDensity tag = endpoint_density(E);
		if (endpoint_type(E) == LEFT) {
			return add_tag(Tags, tag, old_max);
		}
		else {
			return remove_tag(Tags, tag, old_max);
		}
	}

	IntervalList resolve_to_disjoint(IntervalList const& IL) {
		EndpointList EL = convert_to_endpoints(IL);
		IntervalList new_intervals;

		TagHeap tags;
		tags.reserve(IL.size()/2); //this must be enough space

		//Prime the algorithm
		double left = endpoint_location(EL.front());
		IntervalDensity first_tag = endpoint_density(EL.front());
		IntervalDensity old_max = SPARSE;
		tags.push_back(first_tag);

		EndpointList::const_iterator it = std::begin(EL);
		it++;

		for (; it != std::end(EL); it++)
		{
			//If the endpoint results in a new interval
			// in the output stream, then add an interval
			// to the output
			if (AcceptEndpoint(tags, *it, old_max)) {
				//Only add intervals with positive length
				if (endpoint_location(*it) > left || 
					(old_max == DISCONTINUITY && endpoint_density(*it) == DISCONTINUITY)) {
					new_intervals.push_back(TaggedInterval(left, endpoint_location(*it), old_max));
					left = endpoint_location(*it);
				}
			}

		}
		return new_intervals;
	}

	double length(const TaggedInterval &interval)
	{
		return right(interval) - left(interval);
	}

	double total_length(IntervalList const& IL)
	{
		double L = std::accumulate(IL.begin(), IL.end(), 0.0, [](double a, TaggedInterval const& ival) { 
			return a + length(ival); 
		});

		return L;
	}

	IntervalDensity max_density(IntervalDensity d1, IntervalDensity d2)
	{
		return (d1 > d2) ? d1 : d2;
	}

	TaggedInterval intersection(TaggedInterval const& ival1, TaggedInterval const& ival2)
	{
		using std::max; 
		using std::min;
		if (left(ival1) >= right(ival2) || left(ival2) >= right(ival1))
			return TaggedInterval(0, 0, SPARSE);
		else
		{
			return TaggedInterval(max(left(ival1), left(ival2)), 
								  min(right(ival1), right(ival2)), 
								  max_density(density(ival1), density(ival2)));
		}
	}

	bool is_empty(TaggedInterval const& ival)
	{
		if (left(ival) == right(ival))
			return true;
		else
			return false;
	}

	IntervalList	intersect(IntervalList const& IL, TInterval::TaggedInterval const& ival)
	{
		IntervalList	ilist;
		std::transform(IL.begin(), IL.end(), std::back_inserter(ilist), 
			std::bind(intersection, std::placeholders::_1, ival));

		ilist.erase(std::remove_if(ilist.begin(), ilist.end(), is_empty), ilist.end());
		return ilist;
	}
}

