#pragma once
#include "../any.cc"
#include "./any.ii"// testing...






XTAL_ENV_(push)
namespace xtal::math::jonquiere::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

using namespace xtal::math::__test;


////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("any")
{
	TRY_("task")
	{
		TRUE_(true);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
