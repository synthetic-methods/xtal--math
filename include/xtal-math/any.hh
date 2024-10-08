#pragma once
#include <xtal/any.hh>
#include <xtal/all.hh>
#include <Eigen/Dense>

#include "./etc.hh"


XTAL_ENV_(push)
namespace xtal::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class T>
XTAL_REQ duple_q = bond::pack_q<T> and bond::pack_size_n<T> == 2;

template <class U>
//\
XTAL_USE duple_t = algebra::scalar_t<U[2]>;
XTAL_USE duple_t = bond::couple_t<U, U>;


////////////////////////////////////////////////////////////////////////////////

template <int N_arity=-1> requires (N_arity == -1)
XTAL_DEF_(return,inline)
XTAL_LET duple_f(auto &&...xs)
XTAL_0EX -> decltype(auto)
{
	//\
	return algebra::scalar_f(XTAL_REF_(xs)...);
	return bond::couple_f(XTAL_REF_(xs)...);
}
template <int N_arity=-1> requires (N_arity ==  0)
XTAL_DEF_(return,inline)
XTAL_LET duple_f(auto &&x0, auto &&x1, auto &&...xs)
XTAL_0EX -> decltype(auto)
{
	return XTAL_REF_(x0);
}
template <int N_arity=-1> requires (N_arity ==  1)
XTAL_DEF_(return,inline)
XTAL_LET duple_f(auto &&x0, auto &&x1, auto &&...xs)
XTAL_0EX -> decltype(auto)
{
	return duple_f(XTAL_REF_(x0));
}
template <int N_arity=-1> requires (N_arity ==  2)
XTAL_DEF_(return,inline)
XTAL_LET duple_f(auto &&x0, auto &&x1, auto &&...xs)
XTAL_0EX -> decltype(auto)
{
	return duple_f(XTAL_REF_(x0), XTAL_REF_(x1));
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
