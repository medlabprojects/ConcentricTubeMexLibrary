/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include "TupleTransform.h"
#include <boost/mpl/int.hpp>

namespace tuple_transform
{
	struct zip_iterator_tag;
	struct zip_view_tag;

	template <typename Seq, int N>
	struct zip_iterator;

	template <typename Seq1, typename Seq2> struct zip_view;
	template <typename Seq1 > struct zip_with_index_view;

	template <typename T> struct tag_of;

	template <typename Seq1, typename Seq2>
	struct tag_of < zip_view<Seq1, Seq2> >
	{
		typedef zip_view_tag	type;
	};

	template <typename Seq1>
	struct tag_of < zip_with_index_view<Seq1> >
	{
		typedef zip_view_tag	type;
	};

	template <typename Seq, int N>
	struct tag_of < zip_iterator<Seq, N> >
	{
		typedef zip_iterator_tag type;
	};

	template <typename Seq1, typename Seq2>
	struct zip_view
	{
		zip_view(const Seq1& seq1_, const Seq2& seq2_) :
			seq1(seq1_),
			seq2(seq2_)
		{
		}

		typedef Seq1 sequence_type_1;
		typedef Seq2 sequence_type_2;

		template <typename Index>
		struct type_at
		{
			typedef std::tuple< const typename result_of::value_at<Seq1, Index::value>::type&, 
								const typename result_of::value_at<Seq2, Index::value>::type& >	type;
		};

		template <typename Index>
      STRONG_INLINE
		typename type_at<Index>::type	get_item() const
		{
			return std::tie(at<Index::value>(seq1), at<Index::value>(seq2));
		}

		const Seq1 seq1;
		const Seq2 seq2;
	};

	template <typename Seq1>
	struct zip_with_index_view
	{
		zip_with_index_view(const Seq1 &seq1_) :
			seq1(seq1_)
		{
		}

		typedef Seq1 sequence_type_1;

		template <typename Index>
		struct type_at
		{
			typedef std::tuple< const typename result_of::value_at<Seq1, Index::value>::type&,
								 	  const int >	type;
		};

		template <typename Index>
      STRONG_INLINE
		typename type_at<Index>::type	get_item() const
		{
			return typename type_at<Index>::type(at<Index::value>(seq1), Index::value);
		}

		const Seq1 seq1;
	};

	template <typename Seq1, typename Seq2>
   STRONG_INLINE
	zip_view<Seq1, Seq2>	zip(const Seq1& seq1, const Seq2& seq2)
	{
		return zip_view<Seq1, Seq2>(seq1, seq2);
	}

	template <typename Seq1>
   STRONG_INLINE
	zip_with_index_view<Seq1> zip_with_index(const Seq1& seq)
	{
		return zip_with_index_view<Seq1>(seq);
	}

	template <typename Seq, int N>
	struct zip_iterator
	{
		typedef Seq	sequence_type;
		typedef boost::mpl::int_<N>	index;

		zip_iterator(Seq& s_) : s(s_) {}

		Seq& s;
	};

	template <typename tag> struct begin_impl;
	template <typename tag> struct end_impl;
	template <typename tag> struct size_of_impl;
	template <typename tag> struct equal_to_impl;

	template <>
	struct size_of_impl < zip_view_tag >
	{
		template <typename Sequence>
		struct apply
		{
			typedef typename size_of< typename Sequence::sequence_type_1 >::type	type;
		};
	};

	template <>
	struct begin_impl < zip_view_tag >
	{
		template <typename Sequence>
		struct apply
		{
			typedef zip_iterator<Sequence, 0>	type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};
	};

	template <>
	struct at_impl < zip_view_tag >
	{
		template <typename Sequence, int N>
		struct apply
		{
			typedef typename Sequence::template type_at<N>::type	type;

         STRONG_INLINE static type call( Sequence& s )
			{
				return deref(zip_iterator<Sequence, N>(s));
			}
		};

	};

	template <>
	struct end_impl < zip_view_tag >
	{
		template <typename Sequence>
		struct apply
		{
			const static int N = tuple_transform::size_of< Sequence >::type::value;
			typedef zip_iterator<Sequence, N> type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};
	};

	template <>
	struct deref_impl < zip_iterator_tag >
	{
		template <typename It>
		struct apply
		{
			typedef typename It::sequence_type::template type_at<typename It::index>::type	type;
			typedef typename It::index	index;
			
         STRONG_INLINE static type call( const It& it )
			{
				return (it.s).template get_item<index>();
			}
		};
	};

	template <>
	struct next_impl < zip_iterator_tag >
	{
		template <typename It>
		struct apply
		{
			typedef typename It::index	index;
			typedef typename It::sequence_type	sequence_type;
			typedef zip_iterator<sequence_type, index::value+1>	type;

         STRONG_INLINE static type	call( const It& it )
			{
				return type(it.s);
			}
		};
	};

	template <>
	struct equal_to_impl < zip_iterator_tag >
	{
		template <typename It1, typename It2>
		struct apply
		{
			typedef boost::mpl::false_	type;
		};

		template <typename It1>
		struct apply < It1, It1 >
		{
			typedef boost::mpl::true_	type;
		};
	};
}
