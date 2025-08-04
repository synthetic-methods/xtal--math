#pragma once
#include "./any.cc"

#include "../atom/pade/wniplex.hh"



#include "./limit.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("limit")
{
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type; using W_sigma = atom::couple_t<U_sigma[2]>;
	using U_delta = typename U_fit::delta_type; using W_delta = atom::couple_t<U_delta[2]>;
	using U_alpha = typename U_fit::alpha_type; using W_alpha = atom::couple_t<U_alpha[2]>;
	using U_aphex = typename U_fit::aphex_type; using W_aphex = atom::couple_t<U_aphex[2]>;

	auto mt19937_f = typename U_fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	auto constexpr N_half = U_fit::ratio_f(1, 2);
	auto constexpr N_zero = U_fit::ratio_f(0, 2);
	auto constexpr N_min = _std::numeric_limits<U_alpha>::min();
	auto constexpr N_max = _std::numeric_limits<U_alpha>::max();
#ifndef __FINITE_MATH_ONLY__
	auto constexpr N_inf = _std::numeric_limits<U_alpha>::inf();
#endif

	TRY_("wniplex evaluation")
	{
		auto const foo = atom::math::pade::wniplex_f(U_aphex{0.125, 64.00});
		auto const bar = limit_f<3>(foo);
		TRUE_(bond::pack_item_f<1, 0>(foo) <             (bond::pack_item_f<1, 0>(bar)));
		TRUE_(bond::pack_item_f<1, 1>(bar) <             (bond::pack_item_f<1, 1>(foo)));
		TRUE_(bond::pack_item_f<1, 0>(bar) == limit_f< 3>(bond::pack_item_f<1, 0>(foo)));
		TRUE_(bond::pack_item_f<1, 1>(bar) == limit_f<-3>(bond::pack_item_f<1, 1>(foo)));
	}
	TRY_("range evaluation")
	{
	//	W_delta constexpr range{-1, 1};
	//	TRUE_(limit_f<range>(0) == 0);
	//	TRUE_(limit_f<range>(2) == 1);
	//	TRUE_(limit_f<atom::couple_f(-1, 1)>(0) == 0);
	//	TRUE_(limit_f<atom::couple_f(-1, 1)>(2) == 1);

		TRUE_(limit_f<_std::pair{-1, 1}>(0) == 0);
		TRUE_(limit_f<_std::pair{-1, 1}>(2) == 1);

	//	TRUE_(limit_f<_std::pair{-1, 1}>(0.F) == 0.F);
	//	TRUE_(limit_f<_std::pair{-1, 1}>(2.F) == 1.F);

	}
	TRY_("finite evaluation")
	{
		TRUE_(limit_f< 2>(N_min) > N_min);
		TRUE_(limit_f<-2>(N_max) < N_max);
#ifndef __FINITE_MATH_ONLY__
		TRUE_(limit_f<-2>(N_inf) < N_inf);
#endif
	};
	TRY_("threshold evaluation")
	{
		TRUE_(N_min < U_fit::minilon_f(2));
		TRUE_(U_fit::maxilon_f(2) < N_max);

		TRUE_(limit_f<+half>(0.25) == 0.5);
		TRUE_(limit_f<-half>(0.75) == 0.5);

		TRUE_(limit_f<[] XTAL_1FN_(to) (U_fit::upsilon_f(2)*U_fit::haplo_1)>(0.5) > 0.5);
		TRUE_(limit_f<[] XTAL_1FN_(to) (U_fit::dnsilon_f(2)*U_fit::haplo_1)>(0.5) < 0.5);

		TRUE_(limit_f<[] XTAL_1FN_(to) (+U_fit::minilon_f(2))>(N_min) > N_min);
		TRUE_(limit_f<[] XTAL_1FN_(to) (-U_fit::maxilon_f(2))>(N_max) < N_max);
#ifndef __FINITE_MATH_ONLY__
		TRUE_(limit_f<[] XTAL_1FN_(to) (-U_fit::maxilon_f(2))>(N_inf) < N_inf);
#endif

		TRUE_(limit_f<[] XTAL_1FN_(to) (U_fit::upsilon_f(2)*U_fit::minilon_f(2))>(N_min) > N_min);
		TRUE_(limit_f<[] XTAL_1FN_(to) (U_fit::dnsilon_f(2)*U_fit::maxilon_f(2))>(N_max) < N_max);
#ifndef __FINITE_MATH_ONLY__
		TRUE_(limit_f<[] XTAL_1FN_(to) (U_fit::dnsilon_f(2)*U_fit::maxilon_f(2))>(N_inf) < N_inf);
#endif

	};
	TRY_("in-place threshold evaluation via resignation")
	{
		U_alpha o;

	//	Inside:
		o = +0.75; TRUE_(limit_t<+half>::edit_f(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(limit_t<+half>::edit_f(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(limit_t<-half>::edit_f(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(limit_t<-half>::edit_f(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.25; TRUE_(limit_t<+half>::edit_f(o) == +1); TRUE_(o == +0.50);
		o = -0.25; TRUE_(limit_t<+half>::edit_f(o) == -1); TRUE_(o == -0.50);
		o = +0.75; TRUE_(limit_t<-half>::edit_f(o) == +1); TRUE_(o == +0.50);
		o = -0.75; TRUE_(limit_t<-half>::edit_f(o) == -1); TRUE_(o == -0.50);

#ifndef __FINITE_MATH_ONLY__
		o = +inf; TRUE_(limit_t<[] XTAL_1FN_(to) (-U_fit::maxilon_f(1))>::edit_f(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(limit_t<[] XTAL_1FN_(to) (-U_fit::maxilon_f(1))>::edit_f(o) == -1); TRUE_(o  > -inf);
		o = +0.0; TRUE_(limit_t<[] XTAL_1FN_(to) (+U_fit::minilon_f(1))>::edit_f(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(limit_t<[] XTAL_1FN_(to) (+U_fit::minilon_f(1))>::edit_f(o) == -1); TRUE_(o <  -0.0);

		o = inf - inf;
		TRUE_(o != o);
		(void) limit_t<[] XTAL_1FN_(to) (-U_fit::maxilon_f(1))>::edit_f(o);
		TRUE_(o == o);
#endif

	};
	TRY_("in-place threshold evaluation via scaling")
	{
		U_alpha o;

	//	Inside:
		o = +0.75; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::edit_f(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::edit_f(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::edit_f(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::edit_f(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.50; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::edit_f(o) == +1); TRUE_(o  > +0.50);
		o = -0.50; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::haplo_1)>::edit_f(o) == -1); TRUE_(o <  -0.50);
		o = +0.50; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::edit_f(o) == +1); TRUE_(o <  +0.50);
		o = -0.50; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::haplo_1)>::edit_f(o) == -1); TRUE_(o  > -0.50);

#ifndef __FINITE_MATH_ONLY__
		o = +inf; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::maxilon_f(1))>::edit_f(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::dnsilon_f(1)*U_fit::maxilon_f(1))>::edit_f(o) == -1); TRUE_(o  > -inf);
#endif
		o = +0.0; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::minilon_f(1))>::edit_f(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(limit_t<[] XTAL_1FN_(to) (U_fit::upsilon_f(1)*U_fit::minilon_f(1))>::edit_f(o) == -1); TRUE_(o <  -0.0);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
