#pragma once
#include "./any.cc"
#include "./tangy.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using _fit = bond::fit<>;

	using U_sigma = typename _fit::sigma_type;
	using U_delta = typename _fit::delta_type;
	using U_alpha = typename _fit::alpha_type;
	using U_aphex = typename _fit::aphex_type;
	static constexpr U_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr U_alpha two =  2;
	static constexpr U_alpha ten = 10;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(check_f<- 1>(tangy_t< 2>::template method_f< 3>(-0.25),  tanh(-0.25*pie)));
		TRUE_(check_f<-14>(tangy_t< 2>::template method_f< 2>(-0.25),  tanh(-0.25*pie)));
		TRUE_(check_f<-29>(tangy_t< 2>::template method_f< 1>(-0.25),  tanh(-0.25*pie)));
		TRUE_(check_f<-47>(tangy_t< 2>::template method_f< 0>(-0.25),  tanh(-0.25*pie)));

		TRUE_(check_f<- 1>(tangy_t< 1>::template method_f< 3>(-0.25),  tan (-0.25*pie)));
		TRUE_(check_f<-17>(tangy_t< 1>::template method_f< 2>(-0.25),  tan (-0.25*pie)));
		TRUE_(check_f<-32>(tangy_t< 1>::template method_f< 1>(-0.25),  tan (-0.25*pie)));
	//	TRUE_(check_f<- 0>(tangy_t< 1>::template method_f< 0>(-0.25),  tan (-0.25*pie)));

		TRUE_(check_f<-11>(tangy_t<-2>::template method_f< 3>(-0.25), atanh(-0.25)/pie));
		TRUE_(check_f<-19>(tangy_t<-2>::template method_f< 2>(-0.25), atanh(-0.25)/pie));
		TRUE_(check_f<-28>(tangy_t<-2>::template method_f< 1>(-0.25), atanh(-0.25)/pie));
		TRUE_(check_f<-42>(tangy_t<-2>::template method_f< 0>(-0.25), atanh(-0.25)/pie));

		TRUE_(check_f<- 9>(tangy_t<-1>::template method_f< 3>(-0.25), atan (-0.25)/pie));
		TRUE_(check_f<-18>(tangy_t<-1>::template method_f< 2>(-0.25), atan (-0.25)/pie));
		TRUE_(check_f<-27>(tangy_t<-1>::template method_f< 1>(-0.25), atan (-0.25)/pie));
		TRUE_(check_f<-42>(tangy_t<-1>::template method_f< 0>(-0.25), atan (-0.25)/pie));

		for (int i{}; i < 0x100; ++i) {
			U_alpha x_up(+i), x_dn(-i);
			auto const Y_up = atan(x_up)/pie;
			auto const y_up = tangy_t<-1>::template method_f< 2>(x_up);
			auto const Y_dn = atan(x_dn)/pie;
			auto const y_dn = tangy_t<-1>::template method_f< 2>(x_dn);
			TRUE_(check_f<-49>(y_up, Y_up));
		}

//		TRUE_(check_f<-5>(0.5, tangy_t< 2>::template method_f<3>(tangy_t<-2>::template method_f<3>(0.5))));
//		TRUE_(check_f<-6>(0.5, tangy_t< 1>::template method_f<3>(tangy_t<-1>::template method_f<3>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>::template method_f<3>(tangy_t< 2>::template method_f<3>(0.5))));
//		TRUE_(check_f<-4>(0.5, tangy_t<-1>::template method_f<3>(tangy_t< 1>::template method_f<3>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangy_t< 2>::template method_f<2>(tangy_t<-2>::template method_f<2>(0.5))));
//		TRUE_(check_f<-3>(0.5, tangy_t< 1>::template method_f<2>(tangy_t<-1>::template method_f<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>::template method_f<2>(tangy_t< 2>::template method_f<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-1>::template method_f<2>(tangy_t< 1>::template method_f<2>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangy_t< 2>::template method_f<1>(tangy_t<-2>::template method_f<1>(0.5))));
//		TRUE_(check_f<-2>(0.5, tangy_t< 1>::template method_f<1>(tangy_t<-1>::template method_f<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>::template method_f<1>(tangy_t< 2>::template method_f<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-1>::template method_f<1>(tangy_t< 1>::template method_f<1>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangy_t< 2>::template method_f<0>(tangy_t<-2>::template method_f<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t< 1>::template method_f<0>(tangy_t<-1>::template method_f<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>::template method_f<0>(tangy_t< 2>::template method_f<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-1>::template method_f<0>(tangy_t< 1>::template method_f<0>(0.5))));

	};
}
TAG_("tangy trials")
{
	using _fit = bond::fit<>;
	using U_sigma  = typename _fit::sigma_type;
	using U_delta  = typename _fit::delta_type;
	using U_alpha  = typename _fit::alpha_type;
	using U_aphex  = typename _fit::aphex_type;
	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	auto const mt19937_y = half*mt19937_f();

	EST_("tangy< 2; -1> (* Tanh *)\n~! (reals)")
	{
		return tangy_t< 2>::template method_f<-1>(mt19937_y);
	};
	EST_("tangy< 2;  3> (* Tanh *)\n~3 (reals)")
	{
		return tangy_t< 2>::template method_f< 3>(mt19937_y);
	};
	EST_("tangy< 2;  2> (* Tanh *)\n~2 (reals)")
	{
		return tangy_t< 2>::template method_f< 2>(mt19937_y);
	};
	EST_("tangy< 2;  1> (* Tanh *)\n~1 (reals)")
	{
		return tangy_t< 2>::template method_f< 1>(mt19937_y);
	};
	EST_("tangy< 2;  0> (* Tanh *)\n~0 (reals)")
	{
		return tangy_t< 2>::template method_f< 0>(mt19937_y);
	};

	EST_("tangy<-1; -1> (* ArcTan *)\n~! (reals)")
	{
		return tangy_t<-1>::template method_f<-1>(mt19937_y);
	};
	EST_("tangy<-1;  3> (* ArcTan *)\n~3 (reals)")
	{
		return tangy_t<-1>::template method_f< 3>(mt19937_y);
	};
	EST_("tangy<-1;  2> (* ArcTan *)\n~2 (reals)")
	{
		return tangy_t<-1>::template method_f< 2>(mt19937_y);
	};
	EST_("tangy<-1;  1> (* ArcTan *)\n~1 (reals)")
	{
		return tangy_t<-1>::template method_f< 1>(mt19937_y);
	};
	EST_("tangy<-1;  0> (* ArcTan *)\n~0 (reals)")
	{
		return tangy_t<-1>::template method_f< 0>(mt19937_y);
	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
