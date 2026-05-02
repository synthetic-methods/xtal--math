#pragma once
#include "./any.cc"





#include "./tangent.hh"
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
		TRUE_(check_f<-15>(tangent_t< 2>{}.template method< 3>(-2.0),  tanh(-2.0)));
		TRUE_(check_f<-27>(tangent_t< 2>{}.template method< 2>(-2.0),  tanh(-2.0)));
		TRUE_(check_f<-37>(tangent_t< 2>{}.template method< 1>(-2.0),  tanh(-2.0)));
		TRUE_(check_f<-49>(tangent_t< 2>{}.template method< 0>(-2.0),  tanh(-2.0)));

		TRUE_(check_f<- 1>(tangent_t< 2>{}.template method< 3>(-1.0),  tanh(-1.0)));
		TRUE_(check_f<-18>(tangent_t< 2>{}.template method< 2>(-1.0),  tanh(-1.0)));
		TRUE_(check_f<-31>(tangent_t< 2>{}.template method< 1>(-1.0),  tanh(-1.0)));
		TRUE_(check_f<-48>(tangent_t< 2>{}.template method< 0>(-1.0),  tanh(-1.0)));

		TRUE_(check_f<- 1>(tangent_t< 2>{}.template method< 3>(-0.5),  tanh(-0.5)));
		TRUE_(check_f<- 6>(tangent_t< 2>{}.template method< 2>(-0.5),  tanh(-0.5)));
		TRUE_(check_f<-25>(tangent_t< 2>{}.template method< 1>(-0.5),  tanh(-0.5)));
		TRUE_(check_f<-46>(tangent_t< 2>{}.template method< 0>(-0.5),  tanh(-0.5)));

		TRUE_(check_f<- 1>(tangent_t< 1>{}.template method< 3>(-0.5),  tan (-0.5)));
		TRUE_(check_f<- 8>(tangent_t< 1>{}.template method< 2>(-0.5),  tan (-0.5)));
		TRUE_(check_f<-25>(tangent_t< 1>{}.template method< 1>(-0.5),  tan (-0.5)));
		TRUE_(check_f<-47>(tangent_t< 1>{}.template method< 0>(-0.5),  tan (-0.5)));

		TRUE_(check_f<-32>(tangent_t<-2>{}.template method< 3>(-0.5), atanh(-0.5)));
		TRUE_(check_f<-35>(tangent_t<-2>{}.template method< 2>(-0.5), atanh(-0.5)));
		TRUE_(check_f<-40>(tangent_t<-2>{}.template method< 1>(-0.5), atanh(-0.5)));
		TRUE_(check_f<-46>(tangent_t<-2>{}.template method< 0>(-0.5), atanh(-0.5)));

		TRUE_(check_f<-26>(tangent_t<-1>{}.template method< 3>(-0.5), atan (-0.5)));
		TRUE_(check_f<-32>(tangent_t<-1>{}.template method< 2>(-0.5), atan (-0.5)));
		TRUE_(check_f<-37>(tangent_t<-1>{}.template method< 1>(-0.5), atan (-0.5)));
		TRUE_(check_f<-47>(tangent_t<-1>{}.template method< 0>(-0.5), atan (-0.5)));

//		TRUE_(check_f<-5>(0.5, tangent_t< 2>{}.template method<3>(tangent_t<-2>{}.template method<3>(0.5))));
//		TRUE_(check_f<-6>(0.5, tangent_t< 1>{}.template method<3>(tangent_t<-1>{}.template method<3>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>{}.template method<3>(tangent_t< 2>{}.template method<3>(0.5))));
//		TRUE_(check_f<-4>(0.5, tangent_t<-1>{}.template method<3>(tangent_t< 1>{}.template method<3>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangent_t< 2>{}.template method<2>(tangent_t<-2>{}.template method<2>(0.5))));
//		TRUE_(check_f<-3>(0.5, tangent_t< 1>{}.template method<2>(tangent_t<-1>{}.template method<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>{}.template method<2>(tangent_t< 2>{}.template method<2>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-1>{}.template method<2>(tangent_t< 1>{}.template method<2>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangent_t< 2>{}.template method<1>(tangent_t<-2>{}.template method<1>(0.5))));
//		TRUE_(check_f<-2>(0.5, tangent_t< 1>{}.template method<1>(tangent_t<-1>{}.template method<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>{}.template method<1>(tangent_t< 2>{}.template method<1>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-1>{}.template method<1>(tangent_t< 1>{}.template method<1>(0.5))));
//
//		TRUE_(check_f<-1>(0.5, tangent_t< 2>{}.template method<0>(tangent_t<-2>{}.template method<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t< 1>{}.template method<0>(tangent_t<-1>{}.template method<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-2>{}.template method<0>(tangent_t< 2>{}.template method<0>(0.5))));
//		TRUE_(check_f<-1>(0.5, tangent_t<-1>{}.template method<0>(tangent_t< 1>{}.template method<0>(0.5))));

	};
}
TAG_("tangent trials")
{
	using U_fit = bond::fit<>;
	using U_sigma  = typename U_fit::sigma_type;
	using U_delta  = typename U_fit::delta_type;
	using U_alpha  = typename U_fit::alpha_type;
	using U_aphex  = typename U_fit::aphex_type;
	auto mt19937_o = typename U_fit::MT19937{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	auto const mt19937_y = half*mt19937_f();

	EST_("tangent< 2; -1> (* Tanh *)\n~! (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t< 2>{}.template method<-1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent< 2;  3> (* Tanh *)\n~3 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t< 2>{}.template method< 3>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent< 2;  2> (* Tanh *)\n~2 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t< 2>{}.template method< 2>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent< 2;  1> (* Tanh *)\n~1 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t< 2>{}.template method< 1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent< 2;  0> (* Tanh *)\n~0 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t< 2>{}.template method< 0>(half*mt19937_f());
		}
		return w;
	};

	EST_("tangent<-1; -1> (* ArcTan *)\n~! (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t<-1>{}.template method<-1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent<-1;  3> (* ArcTan *)\n~3 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t<-1>{}.template method< 3>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent<-1;  2> (* ArcTan *)\n~2 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t<-1>{}.template method< 2>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent<-1;  1> (* ArcTan *)\n~1 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t<-1>{}.template method< 1>(half*mt19937_f());
		}
		return w;
	};
	EST_("tangent<-1;  0> (* ArcTan *)\n~0 (reals)")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			w *= tangent_t<-1>{}.template method< 0>(half*mt19937_f());
		}
		return w;
	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
