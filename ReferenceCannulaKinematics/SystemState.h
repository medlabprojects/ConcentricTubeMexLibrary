#pragma once

#include "CTRTypedefs.h"
#include "Options.h"
#include "BinaryOperations.h"
#include <boost/mpl/int.hpp>
#include <boost/io/ios_state.hpp>
#include <iomanip>

namespace CTR
{

	template <int N, typename Options>
	struct State;

	namespace internal
	{
		template <typename T>
		struct traits;

		template <int N_, typename Options_>
		struct traits < State<N_, Options_> >
		{
			typedef boost::mpl::int_<N_>		N;
			typedef Options_ Options;
			typedef typename Option::has_option<Options, Option::ComputeGeometry >::type	ComputeGeometry;
			typedef typename Option::has_option<Options, Option::ComputeZta >::type			ComputeZta;
			typedef typename Option::has_option<Options, Option::ExternalLoads >::type		ExternalLoads;
			typedef typename Option::has_option<Options, Option::ComputeZga >::type			ComputeZga;
			typedef typename Option::has_option<Options, Option::ComputeZgf >::type			ComputeZgf;
			typedef typename Option::has_option<Options, Option::ComputeZlf >::type			ComputeZlf;
			typedef typename Option::has_option<Options, Option::ComputeZtf >::type			ComputeZtf;
			typedef typename Option::has_option<Options, Option::ComputeZla >::type			ComputeZla;

			typedef typename boost::mpl::or_< ComputeZta, ComputeZtf >::type				TJac;
			typedef typename boost::mpl::or_< ComputeZla, ComputeZlf >::type				LJac;
			typedef typename boost::mpl::or_< ComputeZgf, ComputeZga >::type				GJac;

			typedef typename boost::mpl::or_< TJac, LJac, GJac >::type						AnyJac;
		};
	}



	template <int N>
	struct ALIGN_SPEC GeometricStates
	{
		typename Vector<3>::type	p;
		typename Vector<4>::type	q;
	};
	struct NoGeometricStates {};

	template <int N>
	struct ALIGN_SPEC LoadingStates
	{
		typename Vector<3>::type	Mt; //!< moment x- and y- components
		typename Vector<3>::type	Nt; //!< internal force 3-vector
	};
	struct NoLoadingStates {};

	template <int N>
	struct ALIGN_SPEC ZtaState
	{
		typename Matrix<2 * N, 2 * N >::type	Zta;
	};
	struct NoZtaState {};

	template <int N>
	struct ALIGN_SPEC ZtfState
	{
		typename Matrix<2 * N, 6>::type		Ztf;
	};
	struct NoZtfState {};

	template <int N>
	struct ALIGN_SPEC ZlaState
	{
		typename Matrix<6, 2 * N >::type	Zla;
	};
	struct NoZlaState {};

	template <int N>
	struct ALIGN_SPEC ZlfState
	{
		typename Matrix<6, 6 >::type		Zlf;
	};
	struct NoZlfState {};

	template <int N>
	struct ALIGN_SPEC ZgfState
	{
		typename Matrix<6, 6 >::type		Zgf;
	};
	struct NoZgfState {};

	template <int N>
	struct ALIGN_SPEC ZgaState
	{
		typename Matrix<6, 2*N >::type		Zga;
	};
	struct NoZgaState {};
	
	/**
		\brief The system state at a single point in space and time
		\param N the number of tubes

		This structure represents the system state at a single
		point in space and time. The state obeys a differential
		state equation, where if \f$y\f$ is a state, we have
		\f[
			\frac{\partial}{\partial s}y = (L_y)_* \xi(s,y)
		\f]

		The entire right hand side of this equation is modeled by 
		CTR::xi_fun
	*/
	template <int N, typename Options>
	struct ALIGN_SPEC State :
		public boost::mpl::if_< typename internal::traits<State<N,Options> >::ComputeGeometry, GeometricStates<N>, NoGeometricStates >::type,
		public boost::mpl::if_< typename internal::traits<State<N,Options> >::ExternalLoads, LoadingStates<N>, NoLoadingStates >::type,
		public boost::mpl::if_< typename internal::traits<State<N,Options> >::ComputeZta, ZtaState<N>, NoZtaState >::type,
		public boost::mpl::if_< typename internal::traits<State<N, Options> >::ComputeZtf, ZtfState<N>, NoZtfState >::type,
		public boost::mpl::if_< typename internal::traits<State<N, Options> >::ComputeZlf, ZlfState<N>, NoZlfState >::type,
		public boost::mpl::if_< typename internal::traits<State<N, Options> >::ComputeZgf, ZgfState<N>, NoZgfState >::type,
		public boost::mpl::if_< typename internal::traits<State<N, Options> >::ComputeZga, ZgaState<N>, NoZgaState >::type,
		public boost::mpl::if_< typename internal::traits<State<N, Options> >::ComputeZla, ZlaState<N>, NoZlaState >::type
	{
		typename Vector<N>::type		Sigma;
		typename Vector<N>::type		Psi; //!< Angles relative to Bishop frame \f$[\psi_{1z},...,\psi_{Nz}]\f$
		typename Vector<N>::type		Mz;  //!< Torsional Moments \f$[m_{1z},...,m_{Nz}]\f$
	};
	
	namespace internal
	{
		const static int print_width = 8;

		template <int N, typename Options, typename Stream>
		void print_impl_t(Stream& stream, const State<N, Options>& state)
		{
			boost::io::ios_all_saver ias(stream);
			using std::endl;
			
			stream << "Sigma = " << std::right << std::fixed << std::setw(print_width) << std::setprecision(5) << state.Sigma.transpose() << std::endl;
			stream << "Psi   = " << std::right << std::fixed << std::setw(print_width) << std::setprecision(5) << state.Psi.transpose() << std::endl;
			stream << "Mz    = " << std::right << std::fixed << std::setw(print_width) << std::setprecision(5) << state.Mz.transpose() << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_g_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream << "p     = " << std::right << std::fixed << std::setw(print_width) << std::setprecision(5) << state.p.transpose() << std::endl;
			stream << "q     = " << std::right << std::fixed << std::setw(print_width) << std::setprecision(5) << state.q.transpose() << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_g_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_g(Stream& stream, const State<N, Options>& state)
		{
			boost::io::ios_all_saver ias(stream);
			print_impl_g_impl(stream, state, typename traits<State<N, Options> >::ComputeGeometry());
		}

		template <int N, typename Options, typename Stream>
		void print_impl_l_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream << "Nt    = " << std::right << std::fixed << std::setw(print_width) << std::setprecision(3) << state.Nt.transpose() << std::endl;
			stream << "Mt   = " << std::right << std::fixed << std::setw(print_width) << std::setprecision(3) << state.Mt.transpose() << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_l_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_l(Stream& stream, const State<N, Options>& state)
		{
			boost::io::ios_all_saver ias(stream);
			print_impl_l_impl(stream, state, typename traits<State<N, Options> >::ExternalLoads());
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zta_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream.precision(3);
			stream << "Zta   = " << std::endl << state.Zta << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zta_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zga_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream.precision(3);
			stream << "Zga   = " << std::endl << state.Zga << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zga_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zlf_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream << "Zlf   = " << std::endl << std::fixed << std::setprecision(3) << state.Zlf << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zlf_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Ztf_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream << "Ztf   = " << std::endl << std::fixed << std::setprecision(3) << state.Ztf << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Ztf_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zgf_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream << "Zgf   = " << std::endl << std::fixed << std::setprecision(3) << state.Zgf << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zgf_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zla_impl(Stream& stream, const State<N, Options>& state, boost::mpl::true_)
		{
			stream << "Zla   = " << std::endl << std::fixed << std::setprecision(3) << state.Zla << std::endl;
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Zla_impl(Stream& stream, const State<N, Options>& state, boost::mpl::false_)
		{
		}

		template <int N, typename Options, typename Stream>
		void print_impl_Z(Stream& stream, const State<N, Options>& state)
		{
			boost::io::ios_all_saver ias(stream);
			print_impl_Zta_impl(stream, state, typename traits<State<N, Options> >::ComputeZta());
			print_impl_Zga_impl(stream, state, typename traits<State<N, Options> >::ComputeZga());
			print_impl_Zgf_impl(stream, state, typename traits<State<N, Options> >::ComputeZgf());
			print_impl_Zla_impl(stream, state, typename traits<State<N, Options> >::ComputeZla());
			print_impl_Ztf_impl(stream, state, typename traits<State<N, Options> >::ComputeZtf());
			print_impl_Zlf_impl(stream, state, typename traits<State<N, Options> >::ComputeZlf());
		}
	}

	template <int N, typename Options, typename Stream>
	Stream& operator<<(Stream& stream, const State<N,Options>& state)
	{
		internal::print_impl_t(stream, state);
		internal::print_impl_g(stream, state);
		internal::print_impl_l(stream, state);
		internal::print_impl_Z(stream, state);

		return stream;
	}

	////////////////////////////////////////////
	//         Binary Operations on state     //
	////////////////////////////////////////////
	namespace internal
	{
		template<int N, typename Options, typename Op>
		void binary_operator_impl_t(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f)
		{
			y.Psi = f(y1.Psi, y2.Psi);
			y.Mz = f(y1.Mz, y2.Mz);
			y.Sigma = f(y1.Sigma, y2.Sigma);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_g(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.p = f(y1.p, y2.p);
			y.q = f(y1.q, y2.q);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_g(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_l(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.Nt = f(y1.Nt, y2.Nt);
			y.Mt = f(y1.Mt, y2.Mt);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_l(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zta(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.Zta = f(y1.Zta, y2.Zta);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zta(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zga(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.Zga = f(y1.Zga, y2.Zga);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zga(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zla(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.Zla = f(y1.Zla, y2.Zla);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zla(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Ztf(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.Ztf = f(y1.Ztf, y2.Ztf);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Ztf(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zgf(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.Zgf = f(y1.Zgf, y2.Zgf);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zgf(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zlf(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::true_)
		{
			y.Zlf = f(y1.Zlf, y2.Zlf);
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Zlf(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f, boost::mpl::false_)
		{
		}

		template<int N, typename Options, typename Op>
		void binary_operator_impl_Z(State<N, Options> &y, const State<N, Options> &y1, const State<N, Options> &y2, const Op& f)
		{
			binary_operator_impl_Zta(y, y1, y2, f, typename Option::has_option<Options, Option::ComputeZta>());
			binary_operator_impl_Zga(y, y1, y2, f, typename Option::has_option<Options, Option::ComputeZga>());
			binary_operator_impl_Zla(y, y1, y2, f, typename Option::has_option<Options, Option::ComputeZla>());
			binary_operator_impl_Ztf(y, y1, y2, f, typename Option::has_option<Options, Option::ComputeZtf>());
			binary_operator_impl_Zgf(y, y1, y2, f, typename Option::has_option<Options, Option::ComputeZgf>());
			binary_operator_impl_Zlf(y, y1, y2, f, typename Option::has_option<Options, Option::ComputeZlf>());
		}

		template<int N, typename Options, typename Op>
		State<N, Options> binary_operator_impl(const State<N, Options>& y1, const State<N, Options>& y2, Op&& f)
		{
			State<N, Options>	y;
			binary_operator_impl_t(y, y1, y2, f);
			binary_operator_impl_g(y, y1, y2, f, typename CTR::Option::has_option<Options, CTR::Option::ComputeGeometry>());
			binary_operator_impl_l(y, y1, y2, f, typename CTR::Option::has_option<Options, CTR::Option::ExternalLoads>());
			binary_operator_impl_Z(y, y1, y2, f);
			return y;
		}
	}

	template <int N, typename Options>
	State<N,Options> operator+(const State<N,Options>& y1, const State<N,Options>& y2)
	{
		return internal::binary_operator_impl(y1, y2, internal::plus());
	}

	template <int N, typename Options>
	State<N, Options> operator-(const State<N, Options>& y1, const State<N, Options>& y2)
	{
		return internal::binary_operator_impl(y1, y2, internal::minus());
	}

	namespace internal
	{
		template<int N, typename Options>
		void scalar_times_impl_t(State<N, Options> &y, double s, const State<N, Options> &y1)
		{
			y.Psi = s*y1.Psi;
			y.Mz = s*y1.Mz;
			y.Sigma = s*y1.Sigma;
#ifdef COUNT_OPS
         g_nmult += N + N + N;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_g(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.p = s*y1.p;
			y.q = s*y1.q;
#ifdef COUNT_OPS
         g_nmult += 7;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_g(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_l(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.Nt = s*y1.Nt;
			y.Mt = s*y1.Mt;
#ifdef COUNT_OPS
         g_nmult += 6;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_l(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_Zta(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.Zta = s*y1.Zta;
#ifdef COUNT_OPS
         g_nmult += 2 * N*N;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_Zta(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_Zga(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.Zga = s*y1.Zga;
#ifdef COUNT_OPS
         g_nmult += 6 * 2 * N;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_Zga(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_Zla(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.Zla = s*y1.Zla;
#ifdef COUNT_OPS
         g_nmult += 2 * N * 6;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_Zla(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_Ztf(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.Ztf = s*y1.Ztf;
#ifdef COUNT_OPS
         g_nmult += 2 * N * 6;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_Ztf(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_Zgf(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.Zgf = s*y1.Zgf;
#ifdef COUNT_OPS
         g_nmult += 6 * 6;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_Zgf(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_Zlf(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::true_)
		{
			y.Zlf = s*y1.Zlf;
#ifdef COUNT_OPS
         g_nmult += 6 * 6;
#endif
		}

		template<int N, typename Options>
		void scalar_times_impl_Zlf(State<N, Options> &y, double s, const State<N,Options> &y1, boost::mpl::false_)
		{
		}

		template<int N, typename Options>
		void scalar_times_impl_Z(State<N, Options> &y, double s, const State<N,Options> &y1)
		{
			scalar_times_impl_Zta(y, s, y1, typename Option::has_option<Options, Option::ComputeZta>());
			scalar_times_impl_Zga(y, s, y1, typename Option::has_option<Options, Option::ComputeZga>());
			scalar_times_impl_Zla(y, s, y1, typename Option::has_option<Options, Option::ComputeZla>());
			scalar_times_impl_Ztf(y, s, y1, typename Option::has_option<Options, Option::ComputeZtf>());
			scalar_times_impl_Zgf(y, s, y1, typename Option::has_option<Options, Option::ComputeZgf>());
			scalar_times_impl_Zlf(y, s, y1, typename Option::has_option<Options, Option::ComputeZlf>());
		}

		template<int N, typename Options>
		State<N, Options> scalar_times_impl(double s, const State<N, Options>& y1)
		{
			State<N, Options>	y;
			scalar_times_impl_t(y, s, y1);
			scalar_times_impl_g(y, s, y1, typename CTR::Option::has_option<Options, CTR::Option::ComputeGeometry>());
			scalar_times_impl_l(y, s, y1, typename CTR::Option::has_option<Options, CTR::Option::ExternalLoads>());
			scalar_times_impl_Z(y, s, y1);
			return y;
		}
	}

	template <int N, typename Options>
	State<N, Options> operator*(double s, const State<N, Options> &y)
	{
		return internal::scalar_times_impl(s, y);
	}

	template <int N, typename Options>
	State<N, Options> operator*(const State<N, Options> &y, double s)
	{
		return internal::scalar_times_impl(s, y);
	}
}
