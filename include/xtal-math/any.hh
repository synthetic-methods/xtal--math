#pragma once
#include <xtal/any.hh>
#include <xtal/all.hh>
#include <Eigen/Dense>

#include "./etc.hh"


XTAL_ENV_(push)
namespace xtal::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//\
template <class U>	using    duple_t = bond::couple_t<U, U>;
template <class U>	using    duple_t = algebra::scalar_t<U[2]>;
template <class T>	concept  duple_q = bond::pack_q<T> and bond::pack_size_n<T> == 2;


////////////////////////////////////////////////////////////////////////////////

template <int N_arity=-1> requires (N_arity == -1)
XTAL_DEF_(short)
XTAL_LET duple_f(auto &&...xs)
noexcept -> decltype(auto)
{
	/*/
	return occur::bundle_f(XTAL_REF_(xs)...);
	/*/
	//\
	return bond::couple_f(XTAL_REF_(xs)...);
	return algebra::scalar_f(XTAL_REF_(xs)...);
	/***/
}
template <int N_arity=-1> requires (N_arity ==  0)
XTAL_DEF_(short)
XTAL_LET duple_f(auto &&x0, auto &&x1, auto &&...xs)
noexcept -> decltype(auto)
{
	return XTAL_REF_(x0);
}
template <int N_arity=-1> requires (N_arity ==  1)
XTAL_DEF_(short)
XTAL_LET duple_f(auto &&x0, auto &&x1, auto &&...xs)
noexcept -> decltype(auto)
{
	return duple_f(XTAL_REF_(x0));
}
template <int N_arity=-1> requires (N_arity ==  2)
XTAL_DEF_(short)
XTAL_LET duple_f(auto &&x0, auto &&x1, auto &&...xs)
noexcept -> decltype(auto)
{
	return duple_f(XTAL_REF_(x0), XTAL_REF_(x1));
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
