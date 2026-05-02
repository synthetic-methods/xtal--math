#pragma once
#include "./any.cc"





#include "./near.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("near")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	auto mt19937_f = typename U_fit::MT19937();
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
 