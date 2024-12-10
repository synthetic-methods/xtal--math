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
	using _op = bond::operating;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;
	XTAL_LET_(T_alpha) one =  1;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;
	XTAL_LET_(T_alpha) inf = _std::numeric_limits<T_alpha>::infinity();

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation")
	{
		TRUE_(cut_f<XTAL_VAL_(0.5F), +1>(0.25) == 0.5);
		TRUE_(cut_f<XTAL_VAL_(0.5F), -1>(0.75) == 0.5);

		TRUE_(cut_f<XTAL_VAL_(+0.5F)>(0.25) == 0.5);
		TRUE_(cut_f<XTAL_VAL_(-0.5F)>(0.75) == 0.5);

		TRUE_(cut_f<XTAL_VAL_(_op::upsilon_1*_op::haplo_1)>(0.5) > 0.5);
		TRUE_(cut_f<XTAL_VAL_(_op::dnsilon_1*_op::haplo_1)>(0.5) < 0.5);

		TRUE_(cut_f<XTAL_VAL_(-_op::maxilon_f(0))>(inf) < inf);
		TRUE_(cut_f<XTAL_VAL_(+_op::minilon_f(0))>(0.0) > 0.0);

		TRUE_(cut_f<XTAL_VAL_(_op::dnsilon_1*_op::maxilon_f(0))>(inf) < inf);
		TRUE_(cut_f<XTAL_VAL_(_op::upsilon_1*_op::minilon_f(0))>(0.0) > 0.0);

	};
	TRY_("in-place evaluation via parameterization")
	{
		T_alpha o;
	
	//	Inside:
		o = +0.75; TRUE_(cut_t<XTAL_VAL_(0.50), +1>::edit(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(cut_t<XTAL_VAL_(0.50), +1>::edit(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(cut_t<XTAL_VAL_(0.50), -1>::edit(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(cut_t<XTAL_VAL_(0.50), -1>::edit(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.25; TRUE_(cut_t<XTAL_VAL_(0.50), +1>::edit(o) == +1); TRUE_(o == +0.50);
		o = -0.25; TRUE_(cut_t<XTAL_VAL_(0.50), +1>::edit(o) == -1); TRUE_(o == -0.50);
		o = +0.75; TRUE_(cut_t<XTAL_VAL_(0.50), -1>::edit(o) == +1); TRUE_(o == +0.50);
		o = -0.75; TRUE_(cut_t<XTAL_VAL_(0.50), -1>::edit(o) == -1); TRUE_(o == -0.50);

	};
	TRY_("in-place evaluation via resignation")
	{
		T_alpha o;

	//	Inside:
		o = +0.75; TRUE_(cut_t<XTAL_VAL_(+0.50)>::edit(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(cut_t<XTAL_VAL_(+0.50)>::edit(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(cut_t<XTAL_VAL_(-0.50)>::edit(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(cut_t<XTAL_VAL_(-0.50)>::edit(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.25; TRUE_(cut_t<XTAL_VAL_(+0.50)>::edit(o) == +1); TRUE_(o == +0.50);
		o = -0.25; TRUE_(cut_t<XTAL_VAL_(+0.50)>::edit(o) == -1); TRUE_(o == -0.50);
		o = +0.75; TRUE_(cut_t<XTAL_VAL_(-0.50)>::edit(o) == +1); TRUE_(o == +0.50);
		o = -0.75; TRUE_(cut_t<XTAL_VAL_(-0.50)>::edit(o) == -1); TRUE_(o == -0.50);

		o = +inf; TRUE_(cut_t<XTAL_VAL_(-_op::maxilon_f(0))>::edit(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(cut_t<XTAL_VAL_(-_op::maxilon_f(0))>::edit(o) == -1); TRUE_(o  > -inf);
		o = +0.0; TRUE_(cut_t<XTAL_VAL_(+_op::minilon_f(0))>::edit(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(cut_t<XTAL_VAL_(+_op::minilon_f(0))>::edit(o) == -1); TRUE_(o <  -0.0);

		o = +inf; TRUE_(cut_t<XTAL_VAL_(-_op::maxilon_f(0))>::edit(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(cut_t<XTAL_VAL_(-_op::maxilon_f(0))>::edit(o) == -1); TRUE_(o  > -inf);
		o = +0.0; TRUE_(cut_t<XTAL_VAL_(+_op::minilon_f(0))>::edit(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(cut_t<XTAL_VAL_(+_op::minilon_f(0))>::edit(o) == -1); TRUE_(o <  -0.0);

	};
	TRY_("in-place evaluation via scaling")
	{
		T_alpha o;

	//	Inside:
		o = +0.75; TRUE_(cut_t<XTAL_VAL_(_op::upsilon_1*_op::haplo_1)>::edit(o) ==  0); TRUE_(o == +0.75);
		o = -0.75; TRUE_(cut_t<XTAL_VAL_(_op::upsilon_1*_op::haplo_1)>::edit(o) ==  0); TRUE_(o == -0.75);
		o = +0.25; TRUE_(cut_t<XTAL_VAL_(_op::dnsilon_1*_op::haplo_1)>::edit(o) ==  0); TRUE_(o == +0.25);
		o = -0.25; TRUE_(cut_t<XTAL_VAL_(_op::dnsilon_1*_op::haplo_1)>::edit(o) ==  0); TRUE_(o == -0.25);

	//	Outside:
		o = +0.50; TRUE_(cut_t<XTAL_VAL_(_op::upsilon_1*_op::haplo_1)>::edit(o) == +1); TRUE_(o  > +0.50);
		o = -0.50; TRUE_(cut_t<XTAL_VAL_(_op::upsilon_1*_op::haplo_1)>::edit(o) == -1); TRUE_(o <  -0.50);
		o = +0.50; TRUE_(cut_t<XTAL_VAL_(_op::dnsilon_1*_op::haplo_1)>::edit(o) == +1); TRUE_(o <  +0.50);
		o = -0.50; TRUE_(cut_t<XTAL_VAL_(_op::dnsilon_1*_op::haplo_1)>::edit(o) == -1); TRUE_(o  > -0.50);

		o = +inf; TRUE_(cut_t<XTAL_VAL_(_op::dnsilon_1*_op::maxilon_f(0))>::edit(o) == +1); TRUE_(o <  +inf);
		o = -inf; TRUE_(cut_t<XTAL_VAL_(_op::dnsilon_1*_op::maxilon_f(0))>::edit(o) == -1); TRUE_(o  > -inf);
		o = +0.0; TRUE_(cut_t<XTAL_VAL_(_op::upsilon_1*_op::minilon_f(0))>::edit(o) == +1); TRUE_(o  > +0.0);
		o = -0.0; TRUE_(cut_t<XTAL_VAL_(_op::upsilon_1*_op::minilon_f(0))>::edit(o) == -1); TRUE_(o <  -0.0);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
