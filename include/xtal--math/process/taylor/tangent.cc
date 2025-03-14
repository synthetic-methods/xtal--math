#pragma once
#include "./any.cc"
#include "./tangent.hh"// testing...

#include "../pade/tangy.hh"



XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;
	static constexpr T_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(check_f<- 1>(tangent_t< 2>::template method_f< 3>(-0.5),  tanh(-0.5)));
		TRUE_(check_f<- 6>(tangent_t< 2>::template method_f< 2>(-0.5),  tanh(-0.5)));
		TRUE_(check_f<-25>(tangent_t< 2>::template method_f< 1>(-0.5),  tanh(-0.5)));
		TRUE_(check_f<-46>(tangent_t< 2>::template method_f< 0>(-0.5),  tanh(-0.5)));

		TRUE_(check_f<- 1>(tangent_t< 1>::template method_f< 3>(-0.5),  tan (-0.5)));
		TRUE_(check_f<- 8>(tangent_t< 1>::template method_f< 2>(-0.5),  tan (-0.5)));
		TRUE_(check_f<-25>(tangent_t< 1>::template method_f< 1>(-0.5),  tan (-0.5)));
		TRUE_(check_f<-47>(tangent_t< 1>::template method_f< 0>(-0.5),  tan (-0.5)));

		TRUE_(check_f<-32>(tangent_t<-2>::template method_f< 3>(-0.5), atanh(-0.5)));
		TRUE_(check_f<-35>(tangent_t<-2>::template method_f< 2>(-0.5), atanh(-0.5)));
		TRUE_(check_f<-40>(tangent_t<-2>::template method_f< 1>(-0.5), atanh(-0.5)));
		TRUE_(check_f<-46>(tangent_t<-2>::template method_f< 0>(-0.5), atanh(-0.5)));

		TRUE_(check_f<-26>(tangent_t<-1>::template method_f< 3>(-0.5), atan (-0.5)));
		TRUE_(check_f<-32>(tangent_t<-1>::template method_f< 2>(-0.5), atan (-0.5)));
		TRUE_(check_f<-37>(tangent_t<-1>::template method_f< 1>(-0.5), atan (-0.5)));
		TRUE_(check_f<-47>(tangent_t<-1>::template method_f< 0>(-0.5), atan (-0.5)));

//		TRUE_(check_f<-5>(0.5, tangent_t< 2>::template method_f<3>(tangent_t<-2>::template method_f<3>(0.5))));
//		TRUE_(check_f<-6>(0.5, tangent_t< 1>::template method_f<3>(tangent_t<-1>::template method_f<3>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>::template method_f<3>(tangent_t< 2>::template method_f<3>(0.5))));
//		TRUE_(check_f<-4>(0.5, tangent_t<-1>::template method_f<3>(tangent_t< 1>::template method_f<3>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangent_t< 2>::template method_f<2>(tangent_t<-2>::template method_f<2>(0.5))));
//		TRUE_(check_f<-3>(0.5, tangent_t< 1>::template method_f<2>(tangent_t<-1>::template method_f<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>::template method_f<2>(tangent_t< 2>::template method_f<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-1>::template method_f<2>(tangent_t< 1>::template method_f<2>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangent_t< 2>::template method_f<1>(tangent_t<-2>::template method_f<1>(0.5))));
//		TRUE_(check_f<-2>(0.5, tangent_t< 1>::template method_f<1>(tangent_t<-1>::template method_f<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>::template method_f<1>(tangent_t< 2>::template method_f<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-1>::template method_f<1>(tangent_t< 1>::template method_f<1>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangent_t< 2>::template method_f<0>(tangent_t<-2>::template method_f<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t< 1>::template method_f<0>(tangent_t<-1>::template method_f<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>::template method_f<0>(tangent_t< 2>::template method_f<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-1>::template method_f<0>(tangent_t< 1>::template method_f<0>(0.5))));

	};
}
TAG_("tangent trials")
{
	using _fit = bond::fit<>;
	using T_sigma  = typename _fit::sigma_type;
	using T_delta  = typename _fit::delta_type;
	using T_alpha  = typename _fit::alpha_type;
	using T_aphex  = typename _fit::aphex_type;
	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	auto const mt19937_y = half*mt19937_f();

	EST_("tangent< 2; -1> (* Tanh *)\n~! (reals)")
	{
		return tangent_t< 2>::template method_f<-1>(mt19937_y);
	};
	EST_("tangent< 2;  3> (* Tanh *)\n~3 (reals)")
	{
		return tangent_t< 2>::template method_f< 3>(mt19937_y);
	};
	EST_("tangent< 2;  2> (* Tanh *)\n~2 (reals)")
	{
		return tangent_t< 2>::template method_f< 2>(mt19937_y);
	};
	EST_("tangent< 2;  1> (* Tanh *)\n~1 (reals)")
	{
		return tangent_t< 2>::template method_f< 1>(mt19937_y);
	};
	EST_("tangent< 2;  0> (* Tanh *)\n~0 (reals)")
	{
		return tangent_t< 2>::template method_f< 0>(mt19937_y);
	};


	EST_("tangent<-1; -1> (* ArcTan *)\n~! (reals)")
	{
		return tangent_t<-1>::template method_f<-1>(mt19937_y);
	};
	EST_("tangent<-1;  3> (* ArcTan *)\n~3 (reals)")
	{
		return tangent_t<-1>::template method_f< 3>(mt19937_y);
	};
	EST_("tangent<-1;  2> (* ArcTan *)\n~2 (reals)")
	{
		return tangent_t<-1>::template method_f< 2>(mt19937_y);
	};
	EST_("tangent<-1;  1> (* ArcTan *)\n~1 (reals)")
	{
		return tangent_t<-1>::template method_f< 1>(mt19937_y);
	};
	EST_("tangent<-1;  0> (* ArcTan *)\n~0 (reals)")
	{
		return tangent_t<-1>::template method_f< 0>(mt19937_y);
	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
