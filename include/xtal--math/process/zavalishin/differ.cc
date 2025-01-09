#pragma once
#include "./any.cc"
#include "./differ.hh"// testing...

#include "./prewarped.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("differ")
{
	using D = differ<>;

	TRY_("instantiation")
	{
		confined_t<D> diff{};
		diff <<= typename D::order_type{1};

		TRUE_(check_f<19>(diff(1), 1));
		TRUE_(check_f<19>(diff(2), 1));
		TRUE_(check_f<19>(diff(4), 2));
		TRUE_(check_f<19>(diff(8), 4));

	}
//	TRY_("scaling")
//	{
//		confined_t<dilate<XTAL_VAL_(bond::operating::patio_2)>, D> diff{};
//		diff <<= typename D::order_type{1};
//
//		TRUE_(check_f<19>(0.15915493667125702, diff(1)));
//		TRUE_(check_f<19>(0.15915493667125702, diff(2)));
//
//	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
