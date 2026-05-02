#pragma once
#include "./any.cc"





#include "./tangy.hh"
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;
	static constexpr U_alpha pie =  3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr U_alpha two =  2;
	static constexpr U_alpha ten = 10;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	auto mt19937_f = typename U_fit::MT19937();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(check_f<- 1>(tangy_t< 2>{}.template method< 3>(-0.25),  tanh(-0.25*pie)));
		TRUE_(check_f<-14>(tangy_t< 2>{}.template method< 2>(-0.25),  tanh(-0.25*pie)));
		TRUE_(check_f<-29>(tangy_t< 2>{}.template method< 1>(-0.25),  tanh(-0.25*pie)));
		TRUE_(check_f<-47>(tangy_t< 2>{}.template method< 0>(-0.25),  tanh(-0.25*pie)));

		TRUE_(check_f<- 1>(tangy_t< 1>{}.template method< 3>(-0.25),  tan (-0.25*pie)));
		TRUE_(check_f<-17>(tangy_t< 1>{}.template method< 2>(-0.25),  tan (-0.25*pie)));
		TRUE_(check_f<-32>(tangy_t< 1>{}.template method< 1>(-0.25),  tan (-0.25*pie)));
	//	TRUE_(check_f<- 0>(tangy_t< 1>{}.template method< 0>(-0.25),  tan (-0.25*pie)));

		TRUE_(check_f<-11>(tangy_t<-2>{}.template method< 3>(-0.25), atanh(-0.25)/pie));
		TRUE_(check_f<-19>(tangy_t<-2>{}.template method< 2>(-0.25), atanh(-0.25)/pie));
		TRUE_(check_f<-28>(tangy_t<-2>{}.template method< 1>(-0.25), atanh(-0.25)/pie));
		TRUE_(check_f<-42>(tangy_t<-2>{}.template method< 0>(-0.25), atanh(-0.25)/pie));

		TRUE_(check_f<- 9>(tangy_t<-1>{}.template method< 3>(-0.25), atan (-0.25)/pie));
		TRUE_(check_f<-18>(tangy_t<-1>{}.template method< 2>(-0.25), atan (-0.25)/pie));
		TRUE_(check_f<-27>(tangy_t<-1>{}.template method< 1>(-0.25), atan (-0.25)/pie));
		TRUE_(check_f<-42>(tangy_t<-1>{}.template method< 0>(-0.25), atan (-0.25)/pie));

		for (int i{}; i < 0x100; ++i) {
			U_alpha x_up(+i), x_dn(-i);
			auto const Y_up = atan(x_up)/pie;
			auto const y_up = tangy_t<-1>{}.template method< 2>(x_up);
			auto const Y_dn = atan(x_dn)/pie;
			auto const y_dn = tangy_t<-1>{}.template method< 2>(x_dn);
			TRUE_(check_f<-49>(y_up, Y_up));
		}

//		TRUE_(check_f<-5>(0.5, tangy_t< 2>{}.template method<3>(tangy_t<-2>{}.template method<3>(0.5))));
//		TRUE_(check_f<-6>(0.5, tangy_t< 1>{}.template method<3>(tangy_t<-1>{}.template method<3>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>{}.template method<3>(tangy_t< 2>{}.template method<3>(0.5))));
//		TRUE_(check_f<-4>(0.5, tangy_t<-1>{}.template method<3>(tangy_t< 1>{}.template method<3>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangy_t< 2>{}.template method<2>(tangy_t<-2>{}.template method<2>(0.5))));
//		TRUE_(check_f<-3>(0.5, tangy_t< 1>{}.template method<2>(tangy_t<-1>{}.template method<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>{}.template method<2>(tangy_t< 2>{}.template method<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-1>{}.template method<2>(tangy_t< 1>{}.template method<2>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangy_t< 2>{}.template method<1>(tangy_t<-2>{}.template method<1>(0.5))));
//		TRUE_(check_f<-2>(0.5, tangy_t< 1>{}.template method<1>(tangy_t<-1>{}.template method<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>{}.template method<1>(tangy_t< 2>{}.template method<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-1>{}.template method<1>(tangy_t< 1>{}.template method<1>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangy_t< 2>{}.template method<0>(tangy_t<-2>{}.template method<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t< 1>{}.template method<0>(tangy_t<-1>{}.template method<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-2>{}.template method<0>(tangy_t< 2>{}.template method<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangy_t<-1>{}.template method<0>(tangy_t< 1>{}.template method<0>(0.5))));

	};
}
TAG_("tangy trials")
{
	using U_fit = bond::fit<>;
	using U_sigma  = typename U_fit::sigma_type;
	using U_delta  = typename U_fit::delta_type;
	using U_alpha  = typename U_fit::alpha_type;
	using U_aphex  = typename U_fit::aphex_type;
	auto mt19937_o = typename U_fit::MT19937{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	EST_("tangy< 2; -1> (* Tanh *)\n~! (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t< 2>{}.template method<-1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy< 2;  3> (* Tanh *)\n~3 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t< 2>{}.template method< 3>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy< 2;  2> (* Tanh *)\n~2 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t< 2>{}.template method< 2>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy< 2;  1> (* Tanh *)\n~1 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t< 2>{}.template method< 1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy< 2;  0> (* Tanh *)\n~0 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t< 2>{}.template method< 0>(half*mt19937_f());
		}
		return w;
	};

	EST_("tangy<-1; -1> (* ArcTan *)\n~! (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t<-1>{}.template method<-1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy<-1;  3> (* ArcTan *)\n~3 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t<-1>{}.template method< 3>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy<-1;  2> (* ArcTan *)\n~2 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t<-1>{}.template method< 2>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy<-1;  1> (* ArcTan *)\n~1 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t<-1>{}.template method< 1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangy<-1;  0> (* ArcTan *)\n~0 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangy_t<-1>{}.template method< 0>(half*mt19937_f());
		}
		return w;
	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
