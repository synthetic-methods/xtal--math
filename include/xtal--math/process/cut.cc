#pragma once
#include "./any.cc"

#include "../atom/pade/uniplex.hh"



#include "./cut.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("cut")
{
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type; using W_sigma = atom::couple_t<U_sigma[2]>;
	using U_delta = typename U_fit::delta_type; using W_delta = atom::couple_t<U_delta[2]>;
	using U_alpha = typename U_fit::alpha_type; using W_alpha = atom::couple_t<U_alpha[2]>;
	using U_aphex = typename U_fit::aphex_type; using W_aphex = atom::couple_t<U_aphex[2]>;

	auto mt19937_f = typename U_fit::MT19937();
	mt19937_f.seed(Catch::rngSeed());

	auto constexpr N_half = U_fit::ratio_f(1, 2);
	auto constexpr N_zero = U_fit::ratio_f(0, 2);
	auto constexpr N_min  =   std::numeric_limits<U_alpha>::min();
	auto constexpr N_max  =   std::numeric_limits<U_alpha>::max();
#ifndef __FINITE_MATH_ONLY__
	auto constexpr N_inf  =   std::numeric_limits<U_alpha>::inf();
#endif

	TRY_("uniplex evaluation")
	{
		auto const foo = atom::math::pade::uniplex_f(U_aphex{0.125, 64.00});
		auto const bar = cut_f<3>(foo);
		TRUE_(bond::pack_item_f<1, 0>(foo) <                (bond::pack_item_f<1, 0>(bar)));
		TRUE_(bond::pack_item_f<1, 1>(bar) <                (bond::pack_item_f<1, 1>(foo)));
		TRUE_(bond::pack_item_f<1, 0>(bar) == cut_f< 3>(bond::pack_item_f<1, 0>(foo)));
		TRUE_(bond::pack_item_f<1, 1>(bar) == cut_f<-3>(bond::pack_item_f<1, 1>(foo)));
	}
	TRY_("range evaluation")
	{
	//	W_delta constexpr range{-1, 1};
	//	TRUE_(cut_f<range>(0) == 0);
	//	TRUE_(cut_f<range>(2) == 1);
	//	TRUE_(cut_f<atom::couple_f(-1, 1)>(0) == 0);
	//	TRUE_(cut_f<atom::couple_f(-1, 1)>(2) == 1);

		TRUE_(cut_f<std::pair{-1, 1}>(0) == 0);
		TRUE_(cut_f<std::pair{-1, 1}>(2) == 1);

	//	TRUE_(cut_f<std::pair{-1, 1}>(0.F) == 0.F);
	//	TRUE_(cut_f<std::pair{-1, 1}>(2.F) == 1.F);

	}
	TRY_("fractional evaluation")
	{
		auto constexpr E_cut = U_fit::ratio_f(1, 3);
		auto constexpr F_cut = cut_f<[] XTAL_1FN_(to) (-E_cut)>;

		TRUE_(check_f<-1>(F_cut( 5.000),  E_cut));
		TRUE_(check_f<-1>(F_cut( 4.000),  E_cut));
		TRUE_(check_f<-1>(F_cut( 3.000),  E_cut));
		TRUE_(check_f<-1>(F_cut( 2.000),  E_cut));
		TRUE_(check_f<-1>(F_cut( 1.000),  E_cut));
		TRUE_(check_f<-1>(F_cut( 0.500),  E_cut));
		TRUE_(check_f<-1>(F_cut( 0.400),  E_cut));
		TRUE_(check_f<-1>(F_cut( 0.300),  0.300));
		TRUE_(check_f<-1>(F_cut( 0.200),  0.200));
		TRUE_(check_f<-1>(F_cut( 0.100),  0.100));
		TRUE_(check_f<-1>(F_cut( 0.000),  0.000));
		TRUE_(check_f<-1>(F_cut(-0.000), -0.000));
		TRUE_(check_f<-1>(F_cut(-0.100), -0.100));
		TRUE_(check_f<-1>(F_cut(-0.200), -0.200));
		TRUE_(check_f<-1>(F_cut(-0.300), -0.300));
		TRUE_(check_f<-1>(F_cut(-0.400), -E_cut));
		TRUE_(check_f<-1>(F_cut(-0.500), -E_cut));
		TRUE_(check_f<-1>(F_cut(-1.000), -E_cut));
		TRUE_(check_f<-1>(F_cut(-2.000), -E_cut));
		TRUE_(check_f<-1>(F_cut(-3.000), -E_cut));
		TRUE_(check_f<-1>(F_cut(-4.000), -E_cut));
		TRUE_(check_f<-1>(F_cut(-5.000), -E_cut));

	};
	TRY_("large evaluation")
	{
		auto constexpr E_cut  = U_fit::diplo_f(0x10);
		auto constexpr F_cut  = cut_f<[] XTAL_1FN_(to) (-E_cut)>;
		auto constexpr F_cut_ = cut_e<[] XTAL_1FN_(to) (-E_cut)>;

		U_alpha x;
		x =  100000.; F_cut_(x); TRUE_(check_f<-1>(x,   E_cut ));
		x =   10000.; F_cut_(x); TRUE_(check_f<-1>(x,   10000.));
		x =    1000.; F_cut_(x); TRUE_(check_f<-1>(x,    1000.));
		x =     100.; F_cut_(x); TRUE_(check_f<-1>(x,     100.));
		x =      10.; F_cut_(x); TRUE_(check_f<-1>(x,      10.));
		x =       1.; F_cut_(x); TRUE_(check_f<-1>(x,       1.));
		x =       0.; F_cut_(x); TRUE_(check_f<-1>(x,       0.));
		x = -     0.; F_cut_(x); TRUE_(check_f<-1>(x, -     0.));
		x = -     1.; F_cut_(x); TRUE_(check_f<-1>(x, -     1.));
		x = -    10.; F_cut_(x); TRUE_(check_f<-1>(x, -    10.));
		x = -   100.; F_cut_(x); TRUE_(check_f<-1>(x, -   100.));
		x = -  1000.; F_cut_(x); TRUE_(check_f<-1>(x, -  1000.));
		x = - 10000.; F_cut_(x); TRUE_(check_f<-1>(x, - 10000.));
		x = -100000.; F_cut_(x); TRUE_(check_f<-1>(x, - E_cut ));

		TRUE_(check_f<-1>(F_cut( 100000.),   E_cut ));
		TRUE_(check_f<-1>(F_cut(  10000.),   10000.));
		TRUE_(check_f<-1>(F_cut(   1000.),    1000.));
		TRUE_(check_f<-1>(F_cut(    100.),     100.));
		TRUE_(check_f<-1>(F_cut(     10.),      10.));
		TRUE_(check_f<-1>(F_cut(      1.),       1.));
		TRUE_(check_f<-1>(F_cut(      0.),       0.));
		TRUE_(check_f<-1>(F_cut(-     0.), -     0.));
		TRUE_(check_f<-1>(F_cut(-     1.), -     1.));
		TRUE_(check_f<-1>(F_cut(-    10.), -    10.));
		TRUE_(check_f<-1>(F_cut(-   100.), -   100.));
		TRUE_(check_f<-1>(F_cut(-  1000.), -  1000.));
		TRUE_(check_f<-1>(F_cut(- 10000.), - 10000.));
		TRUE_(check_f<-1>(F_cut(-100000.), - E_cut ));

	};
	TRY_("finite evaluation")
	{
		TRUE_(cut_f< 2>(N_min) > N_min);
		TRUE_(cut_f<-2>(N_max) < N_max);
#ifndef __FINITE_MATH_ONLY__
		TRUE_(cut_f<-2>(N_inf) < N_inf);
#endif
	};
	TRY_("threshold evaluation")
	{
		TRUE_(N_min < U_fit::minilon_f(2));
		TRUE_(U_fit::maxilon_f(2) < N_max);

		TRUE_(cut_f<+half>(0.25) == 0.5);
		TRUE_(cut_f<-half>(0.75) == 0.5);

		TRUE_(cut_f<[] XTAL_1FN_(to) ( U_fit::upsilon_f(2)*U_fit::haplo_1)>(0.5) > 0.5);
		TRUE_(cut_f<[] XTAL_1FN_(to) ( U_fit::dnsilon_f(2)*U_fit::haplo_1)>(0.5) < 0.5);

		TRUE_(cut_f<[] XTAL_1FN_(to) (+U_fit::minilon_f(2))>(N_min) > N_min);
		TRUE_(cut_f<[] XTAL_1FN_(to) (-U_fit::maxilon_f(2))>(N_max) < N_max);
#ifndef __FINITE_MATH_ONLY__
		TRUE_(cut_f<[] XTAL_1FN_(to) (-U_fit::maxilon_f(2))>(N_inf) < N_inf);
#endif

		TRUE_(cut_f<[] XTAL_1FN_(to) ( U_fit::upsilon_f(2)*U_fit::minilon_f(2))>(N_min) > N_min);
		TRUE_(cut_f<[] XTAL_1FN_(to) ( U_fit::dnsilon_f(2)*U_fit::maxilon_f(2))>(N_max) < N_max);
#ifndef __FINITE_MATH_ONLY__
		TRUE_(cut_f<[] XTAL_1FN_(to) ( U_fit::dnsilon_f(2)*U_fit::maxilon_f(2))>(N_inf) < N_inf);
#endif

	};
	TRY_("in-place threshold evaluation via resignation")
	{
		U_alpha o;

	//	Inside:
		o = +0.75; VOID_(cut_t<+half>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == +0.75);
		o = -0.75; VOID_(cut_t<+half>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == -0.75);
		o = +0.25; VOID_(cut_t<-half>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == +0.25);
		o = -0.25; VOID_(cut_t<-half>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == -0.25);

	//	Outside:
		o = +0.25; VOID_(cut_t<+half>::template method<std::in_place>(o)/* == +1*/); TRUE_(o == +0.50);
		o = -0.25; VOID_(cut_t<+half>::template method<std::in_place>(o)/* == -1*/); TRUE_(o == -0.50);
		o = +0.75; VOID_(cut_t<-half>::template method<std::in_place>(o)/* == +1*/); TRUE_(o == +0.50);
		o = -0.75; VOID_(cut_t<-half>::template method<std::in_place>(o)/* == -1*/); TRUE_(o == -0.50);

#ifndef __FINITE_MATH_ONLY__
		o = +inf; VOID_(cut_t<[] XTAL_1FN_(to) (-U_fit::maxilon_f(1))>::template method<std::in_place>(o)/* == +1*/); TRUE_(o <  +inf);
		o = -inf; VOID_(cut_t<[] XTAL_1FN_(to) (-U_fit::maxilon_f(1))>::template method<std::in_place>(o)/* == -1*/); TRUE_(o  > -inf);
		o = +0.0; VOID_(cut_t<[] XTAL_1FN_(to) (+U_fit::minilon_f(1))>::template method<std::in_place>(o)/* == +1*/); TRUE_(o  > +0.0);
		o = -0.0; VOID_(cut_t<[] XTAL_1FN_(to) (+U_fit::minilon_f(1))>::template method<std::in_place>(o)/* == -1*/); TRUE_(o <  -0.0);

		o = inf - inf;
		TRUE_(o != o);
		(void) cut_t<[] XTAL_1FN_(to) (-U_fit::maxilon_f(1))>::template method<std::in_place>(o);
		TRUE_(o == o);
#endif

	};
	TRY_("in-place threshold evaluation via scaling")
	{
		U_alpha o;

	//	Inside:
		o = +0.75; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == +0.75);
		o = -0.75; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == -0.75);
		o = +0.25; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == +0.25);
		o = -0.25; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* ==  0*/); TRUE_(o == -0.25);

	//	Outside:
		o = +0.50; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* == +1*/); TRUE_(o  > +0.50);
		o = -0.50; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* == -1*/); TRUE_(o <  -0.50);
		o = +0.50; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* == +1*/); TRUE_(o <  +0.50);
		o = -0.50; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::template method<std::in_place>(o)/* == -1*/); TRUE_(o  > -0.50);

#ifndef __FINITE_MATH_ONLY__
		o = +inf; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::maxilon_f(1))>::template method<std::in_place>(o)/* == +1*/); TRUE_(o <  +inf);
		o = -inf; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::maxilon_f(1))>::template method<std::in_place>(o)/* == -1*/); TRUE_(o  > -inf);
#endif
		o = +0.0; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::minilon_f(1))>::template method<std::in_place>(o)/* == +1*/); TRUE_(o  > +0.0);
		o = -0.0; VOID_(cut_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::minilon_f(1))>::template method<std::in_place>(o)/* == -1*/); TRUE_(o <  -0.0);

	};
}
/***/
/**/
TAG_("cut trials")
{
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type; using W_sigma = atom::couple_t<U_sigma[2]>;
	using U_delta = typename U_fit::delta_type; using W_delta = atom::couple_t<U_delta[2]>;
	using U_alpha = typename U_fit::alpha_type; using W_alpha = atom::couple_t<U_alpha[2]>;
	using U_aphex = typename U_fit::aphex_type; using W_aphex = atom::couple_t<U_aphex[2]>;

	auto mt19937_o = typename U_fit::MT19937{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	auto constexpr N_half = U_fit::ratio_f(1, 2);
	auto constexpr N_zero = U_fit::ratio_f(0, 2);
	auto constexpr N_min  =   std::numeric_limits<U_alpha>::min();
	auto constexpr N_max  =   std::numeric_limits<U_alpha>::max();
#ifndef __FINITE_MATH_ONLY__
	auto constexpr N_inf =    std::numeric_limits<U_alpha>::inf();
#endif

	EST_("cut down")
	{
		double w{1};
		for (int i = 0x60; ~--i;) {
			auto x = one/mt19937_f(); (void) cut_t<-8>::template method<std::in_place>(x);
			w += x;
		}
		return w;

	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
