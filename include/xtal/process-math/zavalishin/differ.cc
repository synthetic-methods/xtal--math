#pragma once
#include "./any.cc"
#include "./differ.hh"// testing...

#include "./prewarping.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("differ")
{
	TRY_("instantiation")
	{
		differ_t<> diff{};

		TRUE_(check_f<19>(diff(1), 1));
		TRUE_(check_f<19>(diff(2), 1));
		TRUE_(check_f<19>(diff(4), 2));
		TRUE_(check_f<19>(diff(8), 4));

	}
	TRY_("scaling")
	{
		differ_t<dilate<0, 1>> diff{};

		TRUE_(check_f<19>(0.15915493667125702, diff(1)));
		TRUE_(check_f<19>(0.15915493667125702, diff(2)));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
