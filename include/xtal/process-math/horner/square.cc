#pragma once
#include "./any.cc"
#include "./square.hh"// testing...






XTAL_ENV_(push)
namespace xtal::process::math::horner::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("horner")
{
	TRY_("square")
	{
		TRUE_(square_f(1.0, 2.0, 3.0) == 1.0*1.0 + 2.0*2.0 + 3.0*3.0);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
