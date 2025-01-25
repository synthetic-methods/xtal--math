#pragma once
#include "./any.cc"
#include "./cut.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("cut")
{
	using _fix = bond::fixture<>;
	using T_sigma = typename _fix::sigma_type;
	using T_delta = typename _fix::delta_type;
	using T_alpha = typename _fix::alpha_type;
	using T_aphex = typename _fix::aphex_type;
	T_alpha constexpr inf = _std::numeric_limits<T_alpha>::infinity();

	auto mt19937_f = typename _fix::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation")
	{
		TRUE_(cut_f< half, +1>(0.25) == 0.5);
		TRUE_(cut_f< half, -1>(0.75) == 0.5);

		TRUE_(cut_f<+half>(0.25) == 0.5);
		TRUE_(cut_f<-half>(0.75) == 0.5);

		TRUE_(cut_f<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::haplo_1)>(0.5) > 0.5);
		TRUE_(cut_f<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::haplo_1)>(0.5) < 0.5);

		TRUE_(cut_f<[] XTAL_0FN_(value) (-_fix::maxilon_f(0))>(inf) < inf);
		TRUE_(cut_f<[] XTAL_0FN_(value) (+_fix::minilon_f(0))>(0.0) > 0.0);

		TRUE_(cut_f<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::maxilon_f(0))>(inf) < inf);
		TRUE_(cut_f<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::minilon_f(0))>(0.0) > 0.0);

	};
	TRY_("in-place evaluation via parameterization")
	{
		T_alpha o;
	
	//	Inside:
		o = +0.75; TRUE_(cut_t< half, +1>::static_edit(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(cut_t< half, +1>::static_edit(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(cut_t< half, -1>::static_edit(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(cut_t< half, -1>::static_edit(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.25; TRUE_(cut_t< half, +1>::static_edit(o) == +1); TRUE_(o == +0.50);
		o = -0.25; TRUE_(cut_t< half, +1>::static_edit(o) == -1); TRUE_(o == -0.50);
		o = +0.75; TRUE_(cut_t< half, -1>::static_edit(o) == +1); TRUE_(o == +0.50);
		o = -0.75; TRUE_(cut_t< half, -1>::static_edit(o) == -1); TRUE_(o == -0.50);

	};
	TRY_("in-place evaluation via resignation")
	{
		T_alpha o;

	//	Inside:
		o = +0.75; TRUE_(cut_t<+half>::static_edit(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(cut_t<+half>::static_edit(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(cut_t<-half>::static_edit(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(cut_t<-half>::static_edit(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.25; TRUE_(cut_t<+half>::static_edit(o) == +1); TRUE_(o == +0.50);
		o = -0.25; TRUE_(cut_t<+half>::static_edit(o) == -1); TRUE_(o == -0.50);
		o = +0.75; TRUE_(cut_t<-half>::static_edit(o) == +1); TRUE_(o == +0.50);
		o = -0.75; TRUE_(cut_t<-half>::static_edit(o) == -1); TRUE_(o == -0.50);

		o = +inf; TRUE_(cut_t<[] XTAL_0FN_(value) (-_fix::maxilon_f(0))>::static_edit(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(cut_t<[] XTAL_0FN_(value) (-_fix::maxilon_f(0))>::static_edit(o) == -1); TRUE_(o  > -inf);
		o = +0.0; TRUE_(cut_t<[] XTAL_0FN_(value) (+_fix::minilon_f(0))>::static_edit(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(cut_t<[] XTAL_0FN_(value) (+_fix::minilon_f(0))>::static_edit(o) == -1); TRUE_(o <  -0.0);

		o = +inf; TRUE_(cut_t<[] XTAL_0FN_(value) (-_fix::maxilon_f(0))>::static_edit(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(cut_t<[] XTAL_0FN_(value) (-_fix::maxilon_f(0))>::static_edit(o) == -1); TRUE_(o  > -inf);
		o = +0.0; TRUE_(cut_t<[] XTAL_0FN_(value) (+_fix::minilon_f(0))>::static_edit(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(cut_t<[] XTAL_0FN_(value) (+_fix::minilon_f(0))>::static_edit(o) == -1); TRUE_(o <  -0.0);

	};
	TRY_("in-place evaluation via scaling")
	{
		T_alpha o;

	//	Inside:
		o = +0.75; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::haplo_1)>::static_edit(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::haplo_1)>::static_edit(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::haplo_1)>::static_edit(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::haplo_1)>::static_edit(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.50; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::haplo_1)>::static_edit(o) == +1); TRUE_(o  > +0.50);
		o = -0.50; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::haplo_1)>::static_edit(o) == -1); TRUE_(o <  -0.50);
		o = +0.50; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::haplo_1)>::static_edit(o) == +1); TRUE_(o <  +0.50);
		o = -0.50; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::haplo_1)>::static_edit(o) == -1); TRUE_(o  > -0.50);

		o = +inf; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::maxilon_f(0))>::static_edit(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::dnsilon_1*_fix::maxilon_f(0))>::static_edit(o) == -1); TRUE_(o  > -inf);
		o = +0.0; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::minilon_f(0))>::static_edit(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(cut_t<[] XTAL_0FN_(value) (_fix::upsilon_1*_fix::minilon_f(0))>::static_edit(o) == -1); TRUE_(o <  -0.0);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
