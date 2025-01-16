#pragma once
#include "./any.cc"
#include "./near.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("near")
{
	using T_op = bond::operate<>;
	using T_sigma = typename T_op::sigma_type;
	using T_delta = typename T_op::delta_type;
	using T_alpha = typename T_op::alpha_type;
	using T_aphex = typename T_op::aphex_type;

	auto mt19937_f = typename T_op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("near evaluation")
	{
		TRUE_(near_f(-3.75) == -4.00);
		TRUE_(near_f(-3.50) == -4.00);
		TRUE_(near_f(-3.25) == -4.00);
		TRUE_(near_f(-3.00) == -4.00);
		TRUE_(near_f(-2.75) == -2.00);
		TRUE_(near_f(-2.50) == -2.00);
		TRUE_(near_f(-2.25) == -2.00);
		TRUE_(near_f(-2.00) == -2.00);
		TRUE_(near_f(-1.75) == -2.00);
		TRUE_(near_f(-1.50) == -2.00);
		TRUE_(near_f(-1.25) == -1.00);
		TRUE_(near_f(-1.00) == -1.00);
		TRUE_(near_f(-0.75) == -1.00);
		TRUE_(near_f(-0.50) == -0.50);
		TRUE_(near_f(-0.25) == -0.25);
		TRUE_(near_f(-0.00) == -0.00);
		TRUE_(near_f( 0.00) ==  0.00);
		TRUE_(near_f( 0.25) ==  0.25);
		TRUE_(near_f( 0.50) ==  0.50);
		TRUE_(near_f( 0.75) ==  1.00);
		TRUE_(near_f( 1.00) ==  1.00);
		TRUE_(near_f( 1.25) ==  1.00);
		TRUE_(near_f( 1.50) ==  2.00);
		TRUE_(near_f( 1.75) ==  2.00);
		TRUE_(near_f( 2.00) ==  2.00);
		TRUE_(near_f( 2.25) ==  2.00);
		TRUE_(near_f( 2.50) ==  2.00);
		TRUE_(near_f( 2.75) ==  2.00);
		TRUE_(near_f( 3.00) ==  4.00);
		TRUE_(near_f( 3.25) ==  4.00);
		TRUE_(near_f( 3.50) ==  4.00);
		TRUE_(near_f( 3.75) ==  4.00);
		TRUE_(true);
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
 