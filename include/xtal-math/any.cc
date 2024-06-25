#pragma once
#include "./any.hh"// testing...

#include "./etc.cc"
#include <xtal/etc.hh>




XTAL_ENV_(push)
namespace xtal::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_index>
XTAL_DEF_(return)
XTAL_LET check_f(auto const &u, auto const &v)
XTAL_0EX -> bool
{
	return bond::computrim_f<N_index>(u) == bond::computrim_f<N_index>(v);
}
template <int N_index, int N_limit>
XTAL_DEF_(return)
XTAL_LET check_f(auto const &u, auto const &v)
XTAL_0EX -> int
{
	XTAL_LET Z_index = sign_n<N_index>;
	XTAL_LET Z_limit = sign_n<N_limit>;
	static_assert(Z_index == Z_limit);

	XTAL_IF0
	XTAL_0IF (Z_limit*N_limit <  Z_index*N_index) {
		return check_f<N_limit, N_index>(u, v);
	}
	XTAL_0IF (Z_limit*N_limit == Z_index*N_index) {
		return 0;
	}
	XTAL_0IF (Z_index == -1) {
		return check_f<N_index>(u, v)? N_index: check_f<Z_index + N_index, N_limit>(u, v);
	}
	XTAL_0IF (Z_index ==  1) {
		return check_f<N_limit>(u, v)? N_limit: check_f<N_index, N_limit - Z_limit>(u, v);
	}
}
XTAL_DEF_(return)
XTAL_LET check_f(auto const &u, auto const &v)
XTAL_0EX -> int
{
	return check_f<-1, 1 - (int) bond::operating::fraction.depth>(u, v);
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
