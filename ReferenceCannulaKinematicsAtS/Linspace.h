/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>

namespace Utility
{
	template <typename T>
	class linspace;

	template <typename T>
	class linspace_iterator;

	template <class Iterator> class iterator_traits;
	template <typename T>
	class iterator_traits < linspace_iterator<T> >
	{
	public:
		typedef ptrdiff_t difference_type;
		typedef T	value_type;
		typedef T*	pointer;
		typedef T&	reference;
		typedef std::random_access_iterator_tag	iterator_category;
	};

	template <typename T>
	class linspace_iterator : public std::iterator<std::random_access_iterator_tag, T>
	{
	public:
		typedef linspace_iterator<T>												ThisType;
		typedef typename iterator_traits<ThisType>::difference_type		difference_type;
		typedef typename iterator_traits<ThisType>::value_type			value_type;
		typedef typename iterator_traits<ThisType>::pointer				pointer;
		typedef typename iterator_traits<ThisType>::reference				reference;

		linspace_iterator(const linspace<T>& space, difference_type pos) :
			space_(space),
			pos_(pos)
		{
		}

		linspace_iterator<T>& operator++()
		{
			++pos_;
			return *this;
		}

		linspace_iterator<T> operator++(int)
		{
			linspace_iterator<T> temp = *this;
			++(*this);
			return temp;
		}

		bool operator==(linspace_iterator<T> const& other) const
		{
			//the positions and referenced objects must be the same
			if (this->pos_ == other.pos_ && &(this->space_) == &(other.space_))
			{
				return true;
			}
			else {
				return false;
			}
		}

		bool operator!=(linspace_iterator<T> const& other) const
		{
			return !(this->operator==(other));
		}

		bool operator<(const linspace_iterator<T>& other) const
		{
			return (this->pos_ < other.pos_);
		}

		bool operator>(const linspace_iterator<T>& other) const
		{
			return (this->pos_ > other.pos_);
		}

		bool operator<=(const linspace_iterator<T>& other) const
		{
			return (this->pos_ <= other.pos_);
		}

		bool operator>=(const linspace_iterator<T>& other) const
		{
			return (this->pos_ >= other.pos_);
		}

		linspace_iterator<T>& operator+=(difference_type M)
		{
			pos_ += M;
			return *this;
		}

		linspace_iterator<T>& operator-=(difference_type M)
		{
			pos_ -= M;
			return *this;
		}

		T operator[](difference_type N)
		{
			linspace_iterator<T> i = *this;
			i += N;
			return *i;
		}

		difference_type	operator-(linspace_iterator<T> const& other)
		{
			return pos_ - other.pos_;
		}

		T operator*() const
		{
			double f = static_cast<double>(pos_) / static_cast<double>(space_.GetN() - 1);
			return (1.0 - f)*space_.GetLeft() + f*space_.GetRight();
		}

	private:
		const linspace<T>& space_;
		difference_type pos_;
	};

	template <class T>
	class linspace
	{
	public:
		typedef linspace_iterator<T>	iterator;
		typedef linspace_iterator<T>	const_iterator;

		linspace(T const& left, T const& right, int N) :
			left_(left),
			right_(right),
			N_(N)
		{}

		linspace_iterator<T> begin() const
		{
			return linspace_iterator<T>(*this, 0);
		}

		linspace_iterator<T> end() const
		{
			return linspace_iterator<T>(*this, N_);
		}

		std::reverse_iterator< linspace_iterator<T> >	rbegin() const
		{
			return std::reverse_iterator< linspace_iterator<T> >(end());
		}

		std::reverse_iterator< linspace_iterator<T> >	rend() const
		{
			return std::reverse_iterator< linspace_iterator<T> >(begin());
		}

		const int GetN() const
		{
			return N_;
		}

		const T GetLeft() const
		{
			return left_;
		}

		const T GetRight() const
		{
			return right_;
		}

		const T operator[](ptrdiff_t M)
		{
			return *linspace_iterator<T>(*this, M);
		}

		size_t size() const
		{
			return N_;
		}

		std::vector<T>	toVector() const
		{
			std::vector<T>	v;
			v.reserve(size());

			std::copy(begin(), end(), std::back_inserter(v));
			return v;
		}

	private:
		const T left_;
		const T right_;
		const int N_;
	};
}
