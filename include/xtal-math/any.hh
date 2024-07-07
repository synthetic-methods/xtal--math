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


////////////////////////////////////////////////////////////////////////////////

XTAL_DEF_(return,inline)
XTAL_LET duple_f(auto &&...xs)
XTAL_0EX -> decltype(auto)
{
	//\
	return bond::couple_f(XTAL_REF_(xs)...);
	return algebra::scalar_f(XTAL_REF_(xs)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
