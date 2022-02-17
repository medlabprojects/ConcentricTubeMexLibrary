/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include "TupleTransform.h"

namespace tuple_transform
{
	template <typename SrcSeq>
	struct accumulate_impl
	{
		typedef typename result_of::end<SrcSeq>::type	SrcEnd;

		template <typename I1, typename F, typename T>
      STRONG_INLINE
		static void call(I1 src, const F &f, T &t, mpl::true_) //end
		{}

		template <typename I1, typename F, typename T>
      STRONG_INLINE
		static void call(I1 src, const F &f, T &t, mpl::false_)
		{
			t = f(t, deref(src)); //accumulate in t
			call(tuple_transform::next(src), f, t);
		}

		template <typename I1, typename F, typename T>
      STRONG_INLINE
		static void call(I1 src, const F &f, T &t)
		{
			typename result_of::equal_to<I1, SrcEnd>::type eq;
			return call(src, f, t, eq);
		}

	};

	template <typename Seq, typename F, typename T>
   STRONG_INLINE
	T accumulate(Seq&& s, const F& f, const T& t_)
	{
		T t = t_;
		typedef typename std::remove_reference<Seq>::type	Seq_;
		accumulate_impl<Seq_>::template call(begin(s), f, t);
		return t;
	}

	template <typename SrcSeq>
	struct accumulate_in_place_impl
	{
		typedef typename result_of::end<SrcSeq>::type	SrcEnd;

		template <typename I1, typename F, typename T>
      STRONG_INLINE static void call( I1 src, const F &f, T &t, mpl::true_ ) //end
		{}

		template <typename I1, typename F, typename T>
      STRONG_INLINE static void call( I1 src, const F &f, T &t, mpl::false_ )
		{
			f(t, deref(src)); //accumulate in t (in-place, t is reference type)
			call(tuple_transform::next(src), f, t);
		}

		template <typename I1, typename F, typename T>
      STRONG_INLINE static void call( I1 src, const F &f, T &t )
		{
			typename result_of::equal_to<I1, SrcEnd>::type eq;
			return call(src, f, t, eq);
		}

	};

	template <typename Seq, typename F, typename T>
   STRONG_INLINE
	T& accumulate_in_place(Seq&& s, const F& f, T& t)
	{
		typedef typename std::remove_reference<Seq>::type	Seq_;
		accumulate_in_place_impl<Seq_>::template call(begin(s), f, t);
		return t;
	}
}

