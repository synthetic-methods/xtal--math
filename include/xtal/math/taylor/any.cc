#pragma once
#include "../any.cc"
#include "./any.ii"// testing...






XTAL_ENV_(push)
namespace xtal::math::taylor::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

using namespace xtal::math::__test;


////////////////////////////////////////////////////////////////////////////////

template <int N_trim=0>
XTAL_FN2 check_f(auto const &u, auto const &v)
XTAL_0EX
{
	return bond::computrim_f<N_trim>(u) == bond::computrim_f<N_trim>(v);
}


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
