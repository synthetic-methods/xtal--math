#pragma once
#include "./any.cc"
#include "./termial.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("math")
{
	TRY_("termial")
	{
		TRUE_(termial_f(2.0, 1.0, 2.0, 3.0) == 1.0*1.0 + 2.0*2.0 + 3.0*4.0);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
