/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/bitor.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/not_equal_to.hpp>
#include <boost/mpl/bitand.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/transform.hpp>

namespace CTR
{
	typedef boost::mpl::true_	true_type;
	typedef boost::mpl::false_	false_type;

	namespace Option {
		typedef boost::mpl::int_ < 0x0001 > 	ComputeGeometry;
		static const char* ComputeGeometryName = "Compute Geometry";

		typedef boost::mpl::int_ < 0x0002 >		ExternalLoads;
		static const char* ExternalLoadsName = "Use External Loading";

		typedef boost::mpl::int_ < 0x0004 >		ComputeZta;
		static const char* ComputeZtaName = "Compute Zta matrix";
		
		typedef boost::mpl::int_ < 0x0008 >		ComputeZga;
		static const char* ComputeZgaName = "Compute Zga matrix";

		typedef boost::mpl::int_ < 0x0010 >		ComputeZgf;
		static const char* ComputeZgfName = "Compute Zgf matrix";

		typedef boost::mpl::int_ < 0x0020 >     ComputeZla;
		static const char* ComputeZlaName = "Compute Zla matrix";

		typedef boost::mpl::int_ < 0x0040 >		ComputeZtf;
		static const char* ComputeZtfName = "Compute Ztf matrix";

		typedef boost::mpl::int_ < 0x0080 >		ComputeZlf;
		static const char* ComputeZlfName = "Compute Zlf matrix";

		typedef boost::mpl::int_ < 0x0100 >		ComputeZlaConditional;
		static const char* ComputeZlaConditionalName = "Compute Zla if External Loads are Present";

		typedef boost::mpl::int_ < 0x0000 >		NoOption;
		static const char* NoOptionName = "No additional options";


		//Option Aliases 
		typedef ComputeZta	ComputeStability;
		typedef ComputeZga	ComputeJacobian;
		typedef ComputeZgf	ComputeCompliance;

		/*
			Encodes the dependency graph for the options
		*/
		struct dependent_options
		{
			template <typename Option>
			struct apply
			{
				typedef NoOption type;
			};

		};


			template <>
			struct dependent_options::apply < ComputeGeometry >
			{
				typedef ComputeGeometry type;
			};

			template <>
			struct dependent_options::apply < ExternalLoads >
			{
				typedef boost::mpl::bitor_< ExternalLoads, ComputeGeometry >::type type;
			};

			template <>
			struct dependent_options::apply < ComputeZta >
			{
				typedef boost::mpl::bitor_ < ComputeZta,
											        ComputeZlaConditional > ::type	type;
			};

			template <>
			struct dependent_options::apply < ComputeZla >
			{
				typedef boost::mpl::bitor_ < ComputeZla,
											      ComputeZta, 
											      ExternalLoads,
                                       ComputeZga >::type 	type;
			};

			template <>
			struct dependent_options::apply < ComputeZtf >
			{
				typedef boost::mpl::bitor_ < ComputeZtf,
											 ComputeZlf,
											 ExternalLoads,
                                  ComputeZgf,
                                  ComputeGeometry >::type	type;
			};

			template <>
			struct dependent_options::apply < ComputeZlf >
			{
				typedef boost::mpl::bitor_ < ComputeZlf,
											       ComputeZgf,
                                        ComputeZtf,
											       ExternalLoads,
                                        ComputeGeometry > ::type	type;
			};

			template <>
			struct dependent_options::apply < ComputeZgf >
			{
				typedef boost::mpl::bitor_ < ComputeZgf,
					ComputeGeometry,
					typename apply< ComputeZlf >::type
				> ::type type;
			};

			template <>
			struct dependent_options::apply < ComputeZga >
         {
            typedef boost::mpl::bitor_ < ComputeZga,
                                          ComputeZta,
                                          ComputeZlaConditional,
                                          ComputeGeometry > ::type	type;
			};

		template <typename OptionList, typename Option>
		struct has_option :
			public boost::mpl::if_ <
			typename boost::mpl::not_equal_to <
			typename boost::mpl::bitand_< Option, OptionList >::type,
			boost::mpl::int_<0> > ,
			boost::mpl::true_,
			boost::mpl::false_ >::type
		{};

		template <typename RO>
		struct resolve_conditionals
		{
			typedef typename boost::mpl::if_ <
				boost::mpl::and_< has_option<RO, ComputeZlaConditional>, has_option<RO, ExternalLoads> >,
				typename boost::mpl::bitor_<RO, ComputeZla, ComputeZga >::type,
				RO
				>::type		type;
		};



		template <typename Options, typename Stream>
		void put_options(Stream& S)
		{
			if (has_option<Options, ComputeGeometry>::value)
			{
				S << ComputeGeometryName << std::endl;
			}

			if (has_option<Options, ComputeZta>::value)
			{
				S << ComputeZtaName << std::endl;
			}

			if (has_option<Options, ExternalLoads>::value)
			{
				S << ExternalLoadsName << std::endl;
			}

			if (has_option<Options, ComputeZga>::value)
			{
				S << ComputeZgaName << std::endl;
			}

			if (has_option<Options, ComputeZgf>::value)
			{
				S << ComputeZgfName << std::endl;
			}

			if (has_option<Options, ComputeZlf>::value)
			{
				S << ComputeZlfName << std::endl;
			}

			if (has_option<Options, ComputeZlaConditional>::value)
			{
				S << ComputeZlaConditionalName << std::endl;
			}

			if (has_option<Options, ComputeZla>::value)
			{
				S << ComputeZlaName << std::endl;
			}

			if (has_option<Options, ComputeZtf>::value)
			{
				S << ComputeZtfName << std::endl;
			}
		}
	}

	template <typename O1,
		typename O2 = Option::NoOption,
		typename O3 = Option::NoOption,
		typename O4 = Option::NoOption,
		typename O5 = Option::NoOption,
		typename O6 = Option::NoOption,
		typename O7 = Option::NoOption,
		typename O8 = Option::NoOption>
	struct DeclareOptions
	{
		typedef typename boost::mpl::apply< Option::dependent_options, O1 >::type	O1r;
		typedef typename boost::mpl::apply< Option::dependent_options, O2 >::type	O2r;
		typedef typename boost::mpl::apply< Option::dependent_options, O3 >::type	O3r;
		typedef typename boost::mpl::apply< Option::dependent_options, O4 >::type	O4r;
		typedef typename boost::mpl::apply< Option::dependent_options, O5 >::type	O5r;
		typedef typename boost::mpl::apply< Option::dependent_options, O6 >::type	O6r;
		typedef typename boost::mpl::apply< Option::dependent_options, O7 >::type	O7r;
		typedef typename boost::mpl::apply< Option::dependent_options, O8 >::type	O8r;

		typedef typename boost::mpl::bitor_< O1r, O2r, O3r, O4r, O5r >::type	resolved_options_1_5;
		typedef typename boost::mpl::bitor_< O6r, O7r, O8r >::type				resolved_options_6_8;
		typedef typename boost::mpl::bitor_< resolved_options_1_5, resolved_options_6_8 >::type	resolved_options;

		typedef typename Option::resolve_conditionals<resolved_options>::type	options;
	};

}
