#pragma once
#include "../any.cc"
#include "./any.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::gudermannian::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

using namespace xtal::process::math::_test;


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
