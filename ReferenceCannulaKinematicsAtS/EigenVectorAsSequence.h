/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <Eigen/Dense>
#include <boost/fusion/support/tag_of_fwd.hpp>
#include <boost/fusion/include/tag_of_fwd.hpp>
#include <boost/fusion/iterator.hpp>
#include <boost/mpl/if.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/container/vector/detail/convert_impl.hpp>

template<typename Struct, int Pos>
struct EigenMatrix_iterator
	: boost::fusion::iterator_base<EigenMatrix_iterator<Struct, Pos> >
{
	typedef Struct struct_type;
	typedef boost::mpl::int_<Pos> index;
	typedef boost::fusion::random_access_traversal_tag category;

	EigenMatrix_iterator(Struct& str)
		: struct_(str) {}

	Struct& struct_;
};

namespace boost
{
	namespace fusion
	{
		struct EigenMatrix_tag;
		struct EigenMatrix_iterator_tag;

		namespace traits
		{
			template< typename Scalar, int R, int C, int X, int Y, int Z>
			struct tag_of< Eigen::Matrix<Scalar, R, C, X, Y, Z> >
			{
				typedef EigenMatrix_tag type;
			};

			template< typename Scalar, int R, int C, int X, int Y, int Z>
			struct tag_of< const Eigen::Matrix<Scalar, R, C, X, Y, Z> >
			{
				typedef EigenMatrix_tag type;
			};

			template< typename Struct, int Pos>
			struct tag_of< EigenMatrix_iterator<Struct, Pos> >
			{
				typedef EigenMatrix_iterator_tag type;
			};
		}
	}
}

#include <boost/mpl/bool.hpp>

namespace boost {
	namespace fusion {
		namespace extension
		{
			template<typename>
			struct is_view_impl;

			template<>
			struct is_view_impl<EigenMatrix_tag>
			{
				template<typename Seq>
				struct apply
					: mpl::false_
				{};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template<typename>
			struct is_sequence_impl;

			template<>
			struct is_sequence_impl<EigenMatrix_tag>
			{
				template<typename Seq>
				struct apply
					: mpl::true_
				{};
			};
		}
	}
}

namespace boost {
	namespace fusion
	{
		struct random_access_traversal_tag;

		namespace extension
		{
			template<typename>
			struct category_of_impl;

			template<>
			struct category_of_impl<EigenMatrix_tag>
			{
				template<typename Seq>
				struct apply
				{
					typedef random_access_traversal_tag type;
				};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template<typename>
			struct begin_impl;

			template <>
			struct begin_impl<EigenMatrix_tag>
			{
				template <typename Seq>
				struct apply
				{
					typedef
						::EigenMatrix_iterator<Seq,0>
						type;

					static type
						call(Seq& seq)
					{
						return type(seq);
					}
				};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template <typename>
			struct end_impl;

			template <>
			struct end_impl<EigenMatrix_tag>
			{
				template <typename Seq>
				struct apply
				{
					typedef ::EigenMatrix_iterator<Seq, Seq::SizeAtCompileTime>
						type;

					static type
						call(Seq& seq)
					{
						return type(seq);
					}
				};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template<typename>
			struct size_impl;

			template<>
			struct size_impl<EigenMatrix_tag>
			{
				template<typename Seq>
				struct apply
					: boost::integral_constant<int, Seq::SizeAtCompileTime>
				{};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template<typename>
			struct at_impl;

			template<>
			struct at_impl<EigenMatrix_tag>
			{
				template<typename Seq, typename N>
				struct apply
				{
					typedef typename mpl::if_<
						is_const<Seq>, double, double&>::type type;

					static type
						call(Seq& seq)
					{
						return seq(N::value);
					}
				};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template<typename>
			struct value_at_impl;

			template <>
			struct value_at_impl<EigenMatrix_tag>
			{
				template <typename Seq, typename N>
				struct apply
				{
					typedef double type;
				};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template <typename>
			struct deref_impl;

			template <>
			struct deref_impl<EigenMatrix_iterator_tag>
			{
				template<typename Iterator>
				struct apply;

				template <class Struct, int N>
				struct apply< EigenMatrix_iterator<Struct,N> >
				{
					typedef typename mpl::if_<
						is_const<Struct>, double, double&>::type type;

					static type
						call(EigenMatrix_iterator<Struct,N> const& it)
					{
						return it.struct_(N);
					}
				};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension
		{
			template <typename>
			struct value_of_impl;

			template <>
			struct value_of_impl<EigenMatrix_iterator_tag>
			{
				template <typename It>
				struct apply
				{
					typedef double type;
				};
			};
		}
	}
}

namespace boost {
	namespace fusion {
		namespace extension {
			template <typename>
			struct next_impl;

			template <>
			struct next_impl < EigenMatrix_iterator_tag >
			{
				template<typename Iterator>
				struct apply
				{
					typedef typename Iterator::struct_type struct_type;
					typedef typename Iterator::index index;
					typedef EigenMatrix_iterator<struct_type, index::value + 1> type;

					static type
						call(Iterator const& i)
					{
						return type(i.struct_);
					}
				};
			};
		}
	}
}