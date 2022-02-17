/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <type_traits>
#include <tuple>
#include <boost/mpl/int.hpp>
#include <boost/mpl/if.hpp>

#ifndef STRONG_INLINE
#define STRONG_INLINE 
#endif

namespace Eigen{
	template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
	class Matrix;
}

/**
	\brief Transformation library that allows interoperation of std::tuple with other data types

	This is a template metaprogrammed library which allows for the inter-operation of the std::tuple
	data type with other types such as the Eigen c++ linear algebra library.
*/
namespace tuple_transform {
	namespace mpl = boost::mpl;

	/**
		\brief Tag dispatching mechanism

		All of the types which are compatible with the tuple_transform library are given a 
		"tag" which is simply an empty structure. Each type is given a unique tag, which is used
		for overload resolution to select the correct underlying implementation for each type.
	*/
	template <typename T>
	struct tag_of {};

	/**
		\brief The implementation of the begin() function

		This metafunction is defined for each type using its tag.
	*/
	template <typename Tag>
	struct begin_impl;

	template <typename Tag>
	struct end_impl;

	template <typename Tag>
	struct deref_impl;

	template <typename Tag>
	struct next_impl;

	template <typename Tag>
	struct value_at_impl;

	template <typename Tag>
	struct at_impl;

	template <typename Tag>
	struct equal_to_impl;

	template <typename Tag>
	struct size_of_impl;

	struct std_tuple_tag;
	struct tuple_iterator_tag;

	namespace result_of
	{
		template <typename Sequence, int Index>
		struct value_at :
			public mpl::apply< value_at_impl< typename tag_of<Sequence>::type >, Sequence, mpl::int_<Index> >
		{
		};

		template <typename Sequence>
		struct begin : 
			public mpl::apply< begin_impl< typename tag_of<Sequence>::type >, Sequence>
		{
		};

		template <typename Sequence>
		struct end : 
			public mpl::apply< end_impl< typename tag_of<Sequence>::type >, Sequence >
		{
		};

		template <typename It1, typename It2>
		struct equal_to : 
			public mpl::apply< equal_to_impl< typename tag_of<It1>::type >, It1, It2 >
		{
		};

		template <typename Iterator>
		struct deref :
			public mpl::apply< deref_impl< typename tag_of<Iterator>::type >, 
									  Iterator >
		{
		};

		template <typename Iterator>
		struct next :
			public mpl::apply< next_impl< typename tag_of<Iterator>::type >, Iterator >
		{
		};
	}

	template <typename Sequence>
   STRONG_INLINE
	typename result_of::begin<Sequence>::type
		begin(Sequence& s)
	{
		return mpl::apply< begin_impl< typename tag_of<Sequence>::type >, Sequence>::call(s);
	}

	template <typename Sequence>
   STRONG_INLINE
	typename result_of::end<Sequence>::type
		end(Sequence& s)
	{
		return mpl::apply< end_impl< typename tag_of<Sequence>::type >, Sequence>::call(s);
	}


	template <int M, typename Sequence>
   STRONG_INLINE
	typename at_impl< typename tag_of<Sequence>::type >::template apply<Sequence, M>::type
		at(Sequence& s)
	{
		return at_impl< typename tag_of<Sequence>::type >::template apply<Sequence, M>::call(s);
	}

	template <typename Iterator>
   STRONG_INLINE
	typename result_of::deref<Iterator>::type
		deref(Iterator& it)
	{
		return mpl::apply< deref_impl< typename tag_of<Iterator>::type >, Iterator >::call(it);
	}

	template <typename Iterator>
   STRONG_INLINE
	typename result_of::next<Iterator>::type
		next(Iterator& it)
	{
		return next_impl< typename tag_of<Iterator>::type >::template apply<Iterator>::call(it);
	}

	template <typename Sequence>
	struct size_of :
		boost::mpl::apply< size_of_impl < typename tag_of<Sequence>::type >, Sequence >
	{

	};

	template <typename tuple, int Pos>
	struct tuple_iterator
	{
		typedef mpl::int_<Pos>	index;
		typedef tuple				    tuple_type;

		tuple_iterator(tuple& t_) : t_ref(t_) {}

		tuple& t_ref;
	};

	template <typename ...Args>
	struct tag_of < std::tuple<Args...> >
	{
		typedef std_tuple_tag	type;
	};

	template <typename ...Args>
	struct tag_of < const std::tuple<Args...> >
	{
		typedef std_tuple_tag	type;
	};

	template <typename tuple, int Pos>
	struct tag_of < tuple_iterator<tuple, Pos> >
	{
		typedef tuple_iterator_tag	type;
	};

	template <typename tuple, int Pos>
	struct tag_of < const tuple_iterator<tuple, Pos> >
	{
		typedef tuple_iterator_tag	type;
	};

	template <>
	struct begin_impl < std_tuple_tag >
	{
		template <typename Sequence>
		struct apply
		{
			typedef tuple_iterator<Sequence, 0>	type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};

	};

	template <>
	struct end_impl < std_tuple_tag >
	{
		template <typename Sequence>
		struct apply
		{
			const static int N = std::tuple_size<Sequence>::value;
			typedef tuple_iterator<Sequence, N> type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};
	};

	template <>
	struct equal_to_impl < tuple_iterator_tag >
	{
		template <typename It1, typename It2>
		struct apply
		{
			typedef mpl::false_	type;
		};

		template <typename It1>
		struct apply < It1, It1 >
		{
			typedef mpl::true_ type;
		};
	};

	template <>
	struct value_at_impl < std_tuple_tag >
	{
		template <typename Sequence, typename N>
		struct apply
		{
			typedef typename std::tuple_element<N::value, Sequence>::type type;
		};
	};

	template <>
	struct at_impl < std_tuple_tag >
	{
		template <typename Sequence, int N>
		struct apply
		{
			typedef typename mpl::if_< typename boost::is_const<Sequence>::type,
											  typename std::add_const< typename std::add_lvalue_reference< typename std::tuple_element<N, Sequence>::type >::type >::type,
											  typename std::add_lvalue_reference< typename std::tuple_element<N, Sequence>::type >::type >::type type;

         STRONG_INLINE static type call( Sequence& s )
			{
				return std::get<N>(s);
			}
		};
	};

	template <>
	struct deref_impl < tuple_iterator_tag >
	{
		template <typename Iterator>
		struct apply
		{
			typedef typename Iterator::tuple_type	tuple_type;
			typedef typename Iterator::index		index;
			typedef typename std::add_lvalue_reference< 
										typename std::add_const< 
											typename std::tuple_element<index::value, tuple_type>::type
										>::type
									>::type	type;

         STRONG_INLINE static type call( Iterator& it )
			{
				return std::get<index::value>(it.t_ref);
			}
		};
	};

	template <>
	struct next_impl < tuple_iterator_tag >
	{
		template <typename iterator>
		struct apply
		{
			typedef typename iterator::tuple_type	tuple_type;
			typedef typename iterator::index		index;
			typedef tuple_iterator<tuple_type, index::value + 1 >	type;

         STRONG_INLINE static type	call( iterator& i )
			{
				return type(i.t_ref);
			}
		};
	};

	template <>
	struct size_of_impl < std_tuple_tag >
	{
		template <typename Sequence>
		struct apply
		{
			typedef mpl::int_< std::tuple_size<Sequence>::value >	type;
		};
	};

	struct transform_tag;
	struct transform_tag2;
	struct transform_iterator_tag;

	template <typename F, typename Sequence1, typename Sequence2 = mpl::na> struct transform_view;
	template <typename Sequence, int Pos> struct transform_iterator;

	template <typename F, typename Sequence1>
	struct transform_view< F, Sequence1, mpl::na >
	{
		typedef typename boost::remove_reference<Sequence1>::type	sequence_type;
		typedef typename boost::remove_reference<F>::type			functor_type;

		transform_view(Sequence1& s_, F& f_) : s(s_), f(f_) {}

		template <int N>
      STRONG_INLINE
		typename result_of::value_at< transform_view<F, Sequence1>, N>::type
			transform_call() const
		{
			return f(at<N>(s));
		}

		Sequence1 s;
		F f;
	};

	template <typename F, typename Sequence1, typename Sequence2>
	struct transform_view
	{
		typedef typename boost::remove_reference<Sequence1>::type	sequence_type;
		typedef typename boost::remove_reference<Sequence2>::type	sequence_type2;
		typedef typename boost::remove_reference<F>::type			functor_type;

		transform_view(Sequence1& s_, Sequence2& s2_, F& f_) : s(s_), s2(s2_), f(f_) {}

		template <int N>
      STRONG_INLINE
		typename result_of::value_at< transform_view<F,Sequence1,Sequence2>, N>::type
			transform_call() const
		{
			return f(at<N>(s), at<N>(s2));
		}

		Sequence1 s;
		Sequence2 s2;
		F f;
	};


	//If the sequence passed is an lvalue, the reference can be safely stored inside
	// the transform_view, so make the type const Sequence&
	template <class Sequence, typename F>
	transform_view<const F, 
				   typename std::conditional< std::is_lvalue_reference<Sequence>::value,
											  const Sequence, 
											  typename std::decay<Sequence>::type 
				   >::type 
	>	
		transform(Sequence&& s, const F& f)
	{
		typedef transform_view < const F,
			typename std::conditional < std::is_lvalue_reference<Sequence>::value,
			const Sequence,
			typename std::decay<Sequence>::type
			> ::type
		> ret_type;
		//The constructor, with these types, has the signature
		// transform_view(const Sequence&, const F&)
		return ret_type(s, f);
	}

	template <typename Sequence1, typename Sequence2, typename F>
   STRONG_INLINE
	transform_view<const F,
				   typename std::conditional< std::is_lvalue_reference<Sequence1>::value,
											  const Sequence1,
											  typename std::decay<Sequence1>::type
				   >::type,
				   typename std::conditional< std::is_lvalue_reference<Sequence2>::value,
											  const Sequence2,
											  typename std::decay<Sequence2>::type
				   >::type
	>
		transform(Sequence1&& s, Sequence2&& s2, const F& f)
	{
		typedef transform_view < const F,
			typename std::conditional< std::is_lvalue_reference<Sequence1>::value,
			const Sequence1,
			typename std::decay<Sequence1>::type
			>::type,
			typename std::conditional < std::is_lvalue_reference<Sequence2>::value,
			const Sequence2,
			typename std::decay<Sequence2>::type
			> ::type
		> ret_type;
		return ret_type(s, s2, f);
	}


	template <typename Sequence, int Pos>
	struct transform_iterator
	{
		typedef Sequence				sequence_type;
		typedef mpl::int_<Pos>	index;

		transform_iterator(Sequence& s_) : s(s_){}

		Sequence& s;
	};

	template <typename F, typename Sequence>
	struct tag_of < transform_view<F, Sequence> >
	{
		typedef transform_tag	type;
	};

	template <typename F, typename Sequence>
	struct tag_of < const transform_view<F, Sequence> >
	{
		typedef transform_tag	type;
	};

	template <typename F, typename Sequence1, typename Sequence2>
	struct tag_of < transform_view<F, Sequence1, Sequence2> >
	{
		typedef transform_tag2	type;
	};

	template <typename F, typename Sequence1, typename Sequence2>
	struct tag_of < const transform_view<F, Sequence1, Sequence2> >
	{
		typedef transform_tag2	type;
	};

	template <typename Seq, int N>
	struct tag_of < transform_iterator<Seq, N> >
	{
		typedef transform_iterator_tag	type;
	};

	template <typename Seq, int N>
	struct tag_of < const transform_iterator<Seq, N> >
	{
		typedef transform_iterator_tag	type;
	};

	template <>
	struct begin_impl < transform_tag >
	{
		template <typename Sequence>
		struct apply
		{
			typedef transform_iterator<Sequence, 0>	type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};
	};

	template <>
	struct begin_impl < transform_tag2 >
	{
		template <typename Sequence>
		struct apply
		{
			typedef transform_iterator<Sequence, 0>	type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};
	};

	template <>
	struct end_impl < transform_tag >
	{
		template <typename Sequence>
		struct apply
		{
			typedef typename size_of<Sequence>::type	size;
			typedef transform_iterator<Sequence, size::value>	type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};
	};

	template <>
	struct end_impl < transform_tag2 >
	{
		template <typename Sequence>
		struct apply
		{
			typedef typename size_of<Sequence>::type	size;
			typedef transform_iterator<Sequence, size::value>	type;

         STRONG_INLINE static type	call( Sequence& seq )
			{
				return type(seq);
			}
		};
	};

	template <>
	struct next_impl < transform_iterator_tag >
	{
		template <typename Iterator>
		struct apply
		{
			typedef typename Iterator::index	index;
			typedef typename Iterator::sequence_type	sequence_type;
			typedef transform_iterator< sequence_type, index::value + 1>	type;

         STRONG_INLINE static type call( Iterator& it )
			{
				return type(it.s);
			}
		};
	};

	template <>
	struct deref_impl < transform_iterator_tag >
	{
		template <typename Iterator>
		struct apply
		{
			typedef typename Iterator::index				index;
			typedef typename Iterator::sequence_type		sequence_type;

			typedef typename result_of::value_at<sequence_type, index::value >::type	type;

         STRONG_INLINE static type call( Iterator& it )
			{
				return it.s.template transform_call<index::value>();
			}
		};
	};

	template <>
	struct equal_to_impl < transform_iterator_tag >
	{
		template <typename It1, typename It2 >
		struct apply
		{
			typedef mpl::false_ type;
		};

		template <typename It1>
		struct apply < It1, It1 >
		{
			typedef mpl::true_ type;
		};
	};

	template <>
	struct value_at_impl < transform_tag >
	{
		template <typename Sequence, typename N>
		struct apply
		{
			typedef typename Sequence::functor_type													functor_type;
			typedef typename Sequence::sequence_type												   sequence_type;
			typedef typename result_of::value_at<sequence_type, N::value>::type				value_type;
			typedef typename std::result_of< functor_type(const value_type&) >::type		type;
		};
	};

	template <>
	struct value_at_impl < transform_tag2 >
	{
		template <typename Sequence, typename N>
		struct apply
		{
			typedef typename Sequence::functor_type	functor_type;
			typedef typename Sequence::sequence_type	sequence_type;
			typedef typename Sequence::sequence_type2	sequence_type2;
			typedef typename result_of::value_at<sequence_type, N::value>::type		value_type;
			typedef typename result_of::value_at<sequence_type2, N::value>::type	value_type2;
			typedef typename std::result_of< functor_type(const value_type&, const value_type2&) >::type	type;
		};
	};

	template <>
	struct at_impl < transform_tag >
	{
		template <typename Sequence, int N>
		struct apply
		{
			typedef typename Sequence::functor_type	functor_type;
			typedef typename Sequence::sequence_type	sequence_type;
			typedef typename result_of::value_at<sequence_type, N>::type				value_type;
			typedef typename std::result_of< functor_type(const value_type&) >::type	type;

         STRONG_INLINE static type call( Sequence& seq )
			{
				return seq.f(at<N>(seq.s));
			}
		};
	};

	template <>
	struct at_impl < transform_tag2 >
	{
		template <typename Sequence, int N>
		struct apply
		{
			typedef typename Sequence::functor_type	functor_type;
			typedef typename Sequence::sequence_type	sequence_type;
			typedef typename Sequence::sequence_type2	sequence_type2;
			typedef typename result_of::value_at<sequence_type, N>::type		value_type;
			typedef typename result_of::value_at<sequence_type2, N>::type		value_type2;
			typedef typename std::result_of< functor_type(const value_type&, const value_type2&) >::type	type;

         STRONG_INLINE static type call( Sequence& seq )
			{
				return seq.f(at<N>(seq.s), at<N>(seq.s2));
			}
		};
	};

	template <>
	struct size_of_impl < transform_tag >
	{
		template <typename Sequence>
		struct apply
		{
			typedef typename Sequence::sequence_type	sequence_type;
			typedef typename size_of<sequence_type>::type	type;
		};
	};

	template <>
	struct size_of_impl < transform_tag2 >
	{
		template <typename Sequence>
		struct apply
		{
			typedef typename Sequence::sequence_type	sequence_type;
			typedef typename size_of<sequence_type>::type	type;
		};
	};

	template <typename SrcSeq, typename DestSeq>
	struct copy_impl
	{
		typedef typename result_of::end<SrcSeq>::type	SrcEnd;

		template <typename I1, typename I2>
      STRONG_INLINE
		static void call(I1 src, I2 dest, mpl::true_)
		{}

		template <typename I1, typename I2>
      STRONG_INLINE
		static void call(I1 src, I2 dest, mpl::false_)
		{
			tuple_transform::deref(dest) = tuple_transform::deref(src);
			call(tuple_transform::next(src), tuple_transform::next(dest));
		}

		template <typename I1, typename I2>
      STRONG_INLINE
		static void call(I1 src, I2 dest)
		{
			typename result_of::equal_to<I1, SrcEnd>::type eq;
			return call(src, dest, eq);
		}

	};

	template <typename SrcSeq, typename DestSeq>
   STRONG_INLINE
	void copy(const SrcSeq& s, DestSeq&& d)
	{
		typedef typename std::remove_reference<DestSeq>::type	DestType;
		copy_impl<const SrcSeq, DestType>::template call(begin(s), begin(d));
	}

	//Declare the types for the tag names
	struct EigenMatrix_tag;
	struct EigenMatrix_iterator_tag;

	template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
	struct tag_of < Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
	{
		typedef EigenMatrix_tag	type;
	};

	template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
	struct tag_of < const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
	{
		typedef EigenMatrix_tag	type;
	};

	/**
		\brief Allows Eigen matrices to work directly with tuple_transform functions.

		This iterator class allows Eigen matrix types to be passed in as values
		to tuple_transform calls. As an example, consider the case where we may
		have several compound types Tube1, Tube2, and Tube3. We can now call a two-argument
		transformation which takes each Tube variable, pairs it with a double from an Eigen vector
		type, and calls a two argument function.  For example, take the following function object.

		~~~{.cpp}
			struct thing_to_do
			{
				template <typename Tube>
				double operator()(Tube const& t, double d)
				{
					return t.doSomething(d);
				}
			};
		~~~

		This function object can be applied to a tuple of heterogeneous Tube types, as follows:

		~~~{.cpp}
			Tube1 tube_1;
			Tube2 tube_2;
			Tube3 tube_3;
			std::tuple<Tube1, Tube2, Tube3>	tup = std::make_tuple(tube_1, tube_2, tube_2);
			Eigen::Vector3d args; args[0] = 1.0; args[1] = 2.0; args[2] = 3.0;
			auto results = tuple_transform::transform(tup, args, thing_to_do());

			Eigen::Vector3d results_vector;
			tuple_transform::copy( results, results_vector);

			std::cout << results_vector << std::endl;
		~~~

		Note the last copy operation. This is necessary to actually evaluate the transformation, as
		the "results" variable is a lazy transformation object that simply contains instructions for
		how to compute the results. The copy operation forces the results to be computed, after which
		they can be displayed.
	*/
	template <typename _matrix_type, int N >
	struct EigenMatrix_iterator
	{
		typedef _matrix_type							matrix_type;
		typedef typename matrix_type::Scalar	scalar_type;
		typedef mpl::int_<N>							index;

		EigenMatrix_iterator(matrix_type& m_) : m(m_) {}

		matrix_type& m;
	};

	template <typename _matrix_type, int N>
	struct tag_of < EigenMatrix_iterator<_matrix_type, N> >
	{
		typedef EigenMatrix_iterator_tag	type;
	};

	template <typename _matrix_type, int N>
	struct tag_of < const EigenMatrix_iterator<_matrix_type, N> >
	{
		typedef EigenMatrix_iterator_tag	type;
	};

	template <>
	struct next_impl < EigenMatrix_iterator_tag >
	{
		template <typename Iterator>
		struct apply
		{
			typedef typename Iterator::index							index;
			typedef typename Iterator::matrix_type						matrix_type;
			typedef EigenMatrix_iterator<matrix_type, index::value + 1>	type;

         STRONG_INLINE static type call( Iterator it )
			{
				return type(it.m);
			}
		};
	};

	template <>
	struct size_of_impl < EigenMatrix_tag >
	{
		template <typename Matrix>
		struct apply
		{
			typedef mpl::int_<Matrix::RowsAtCompileTime>	rows;
			typedef mpl::int_<Matrix::ColsAtCompileTime>	cols;
			typedef typename mpl::int_< rows::value * cols::value >	type;
		};
	};

	template <>
	struct deref_impl < EigenMatrix_iterator_tag >
	{
		template <typename Iterator>
		struct apply
		{
			typedef typename Iterator::scalar_type	scalar_type;
			typedef typename Iterator::matrix_type	matrix_type;
			typedef typename mpl::if_ < typename boost::is_const<matrix_type>::type,
												const  scalar_type&,
												scalar_type& > ::type	type;
			
         STRONG_INLINE static type	call( Iterator const& it )
			{
				return it.m[Iterator::index::value];
			}
		};
	};

	template <>
	struct equal_to_impl < EigenMatrix_iterator_tag >
	{
		template <typename It1, typename It2>
		struct apply
		{
			typedef mpl::false_ type;
		};

		template <typename It1>
		struct apply < It1, It1 >
		{
			typedef mpl::true_ type;
		};
	};

	template <>
	struct begin_impl < EigenMatrix_tag >
	{
		template <typename Matrix>
		struct apply
		{
			typedef EigenMatrix_iterator< Matrix, 0 >	type;

         STRONG_INLINE static type call( Matrix& m )
			{
				return type(m);
			}
		};
	};

	template <>
	struct end_impl < EigenMatrix_tag >
	{
		template <typename Matrix>
		struct apply
		{
			typedef mpl::int_<Matrix::RowsAtCompileTime>	rows;
			typedef mpl::int_<Matrix::ColsAtCompileTime>	cols;
			typedef typename mpl::int_< rows::value * cols::value >	N;
			typedef EigenMatrix_iterator< Matrix, N::value >	type;
			

         STRONG_INLINE static type call( Matrix& m )
			{
				return type(m);
			}
		};
	};

	template <>
	struct value_at_impl < EigenMatrix_tag >
	{
		template <typename Matrix, typename N>
		struct apply
		{
			typedef typename Matrix::Scalar	scalar_type;
			typedef scalar_type	type;
		};
	};

	template <>
	struct at_impl < EigenMatrix_tag >
	{
		template <typename Matrix, int N>
		struct apply
		{
			typedef typename Matrix::Scalar scalar_type;
			typedef typename mpl::if_< typename boost::is_const< Matrix >::type,
											  const scalar_type&,
											  scalar_type& >::type type;

         STRONG_INLINE static type call( Matrix& m )
			{
				return m(N);
			}
		};
	};

	template <typename view, int N>
	struct EigenMatrix_view_iterator
	{
		EigenMatrix_view_iterator(view& v_) : v(v_) {}

		typedef view						view_type;
		typedef mpl::int_<N>			index;

		view_type& v;
	};

	template <typename matrix>
	struct EigenMatrix_column_view
	{
		typedef mpl::int_< matrix::ColsAtCompileTime >	size;
		typedef matrix	matrix_type;
		typedef decltype(std::declval<matrix>().col(std::declval<int>()))	expr_type;

		EigenMatrix_column_view(matrix& m_) : m(m_) {}

      STRONG_INLINE
		const expr_type
			view(int i) const
		{
			return m.col(i);
		}

      STRONG_INLINE
		expr_type
			view(int i)
		{
			return m.col(i);
		}

		matrix& m;
	};

	template <typename matrix>
	struct EigenMatrix_row_view
	{
		typedef mpl::int_< matrix::RowsAtCompileTime >	size;
		typedef matrix matrix_type;
		typedef decltype(std::declval<matrix>().row(std::declval<int>()))	expr_type;

		EigenMatrix_row_view(matrix& m_) : m(m_) {}

      STRONG_INLINE
		const expr_type
			view(int i) const
		{
			return m.row(i);
		}

      STRONG_INLINE
		expr_type
			view(int i) 
		{
			return m.row(i);
		}

		matrix& m;
	};

	struct EigenMatrix_column_view_tag;
	struct EigenMatrix_row_view_tag;
	struct EigenMatrix_view_iterator_tag;

	template <typename view, int N>
	struct tag_of < EigenMatrix_view_iterator<view, N> >
	{
		typedef EigenMatrix_view_iterator_tag	type;
	};

	template <typename view, int N>
	struct tag_of < const EigenMatrix_view_iterator<view, N> >
	{
		typedef EigenMatrix_view_iterator_tag	type;
	};

	template <typename matrix>
	struct tag_of < EigenMatrix_column_view<matrix> >
	{
		typedef EigenMatrix_column_view_tag	type;
	};

	template <typename matrix>
	struct tag_of < const EigenMatrix_column_view<matrix> >
	{
		typedef EigenMatrix_column_view_tag	type;
	};

	template <typename matrix>
	struct tag_of < EigenMatrix_row_view<matrix> >
	{
		typedef EigenMatrix_row_view_tag	type;
	};

	template <typename matrix>
	struct tag_of < const EigenMatrix_row_view<matrix> >
	{
		typedef EigenMatrix_row_view_tag	type;
	};

	template <>
	struct deref_impl < EigenMatrix_view_iterator_tag >
	{
		template <typename Iterator>
		struct apply
		{
			typedef typename Iterator::view_type	view_type;
			typedef typename view_type::matrix_type	matrix_type;
			typedef typename Iterator::index		index;
			typedef typename view_type::expr_type	type;

         STRONG_INLINE static type call( Iterator& it )
			{
				return it.v.view(index::value);
			}
		};
	};

	template <>
	struct next_impl < EigenMatrix_view_iterator_tag >
	{
		template <typename Iterator>
		struct apply
		{
			typedef typename Iterator::view_type	view_type;
			typedef typename Iterator::index		index;
			typedef EigenMatrix_view_iterator<view_type, index::value + 1>	type;

         STRONG_INLINE static type call( Iterator& it )
			{
				return type(it.v);
			}
		};
	};

	template <>
	struct begin_impl < EigenMatrix_column_view_tag >
	{
		template <typename View>
		struct apply
		{
			typedef EigenMatrix_view_iterator<View, 0>	type;

         STRONG_INLINE static type call( View& v )
			{
				return type(v);
			}
		};
	};

	template <>
	struct end_impl < EigenMatrix_column_view_tag >
	{
		template <typename View>
		struct apply
		{
			typedef typename View::size									N;
			typedef EigenMatrix_view_iterator<View, N::value> type;

         STRONG_INLINE static type call( View& v )
			{
				return type(v);
			}
		};
	};

	template <>
	struct at_impl < EigenMatrix_column_view_tag >
	{
		template <typename View, int N>
		struct apply
		{
			typedef	typename View::expr_type	type;

         STRONG_INLINE static type call( View& v )
			{
				return v.view(N);
			}
		};
	};

	template <>
	struct value_at_impl < EigenMatrix_column_view_tag >
	{
		template <typename View, typename N>
		struct apply
		{
			typedef	typename View::expr_type	type;
		};
	};

	template <>
	struct size_of_impl < EigenMatrix_column_view_tag >
	{
		template <typename View>
		struct apply
		{
			typedef typename View::size	type;
		};
	};

	template <>
	struct begin_impl < EigenMatrix_row_view_tag >
	{
		template <typename View>
		struct apply
		{
			typedef EigenMatrix_view_iterator<View, 0>	type;

         STRONG_INLINE static type call( View& v )
			{
				return type(v);
			}
		};
	};

	template <>
	struct end_impl < EigenMatrix_row_view_tag >
	{
		template <typename View>
		struct apply
		{
			typedef typename View::size									N;
			typedef EigenMatrix_view_iterator<View, N::value> type;

         STRONG_INLINE static type call( View& v )
			{
				return type(v);
			}
		};
	};

	/**
		\brief The implementation of at for row views of Eigen matrices
	*/
	template <>
	struct at_impl < EigenMatrix_row_view_tag >
	{
		template <typename View, int N>
		struct apply
		{
			typedef typename View::expr_type	type;

         STRONG_INLINE static type call( View& v )
			{
				return v.view(N);
			}
		};
	};

	/**
		\brief The implementation of value_at for row views of Eigen matrices
	*/
	template <>
	struct value_at_impl < EigenMatrix_row_view_tag >
	{
		template <typename View, typename N>
		struct apply
		{
			typedef	typename View::expr_type	type;
		};
	};
	
	/**
		\brief The implementation of size_of for row views of Eigen matrices
	*/
	template <>
	struct size_of_impl < EigenMatrix_row_view_tag >
	{
		//Apply implementation just pulls the size
		// from the view itself
		template <typename View>
		struct apply
		{
			typedef typename View::size	type;
		};
	};

	/**
		\brief Return a view of the columns of a matrix

		This function takes an Eigen matrix type and returns a view
		into this matrix by columns, for use with tuple_transform::transform
		and tuple_transform::copy.  For example, as a trivial example, consider

		~~~{.cpp}
			std::tuple<double,double,double>	tup = std::make_tuple(0.0,1.0,2.0);
			Eigen::Matrix<double,1,3>	row_vector;
			tuple_transform::copy( tup, tuple_transform::columns_of(row_vector) );
		~~~

		This code will place the tuple values (0.0,1.0,2.0) into the corresponding
		entries in the row vector. 
	*/
	template <typename matrix>
   STRONG_INLINE
	EigenMatrix_column_view<matrix>	columns_of(matrix &m)
	{
		return EigenMatrix_column_view<matrix>(m);
	}

	/**
		\brief Return a view of the rows of a matrix

		This function takes an Eigen matrix type and returns a view
		into this matrix by columns, for use with tuple_transform::transform
		and tuple_transform::copy.  For example, as a trivial example, consider

		~~~{.cpp}
			std::tuple<double,double,double>	tup = std::make_tuple(0.0,1.0,2.0);
			Eigen::Matrix<double,3,1>	column_vector;
			tuple_transform::copy( tup, tuple_transform::rows_of(column_vector) );
		~~~

		This code will place the tuple values (0.0,1.0,2.0) into the corresponding
		entries in the column vector.
	*/
	template <typename matrix>
   STRONG_INLINE
	EigenMatrix_row_view<matrix> rows_of(matrix &m)
	{
		return EigenMatrix_row_view<matrix>(m);
	}
}
