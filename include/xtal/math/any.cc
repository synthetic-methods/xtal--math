#pragma once
#include "./any.hh"// testing...

#include <Eigen/Dense>
#include <catch2/catch_all.hpp>

#define UNTRUE_(...)   REQUIRE(not (__VA_ARGS__))
#define   TRUE_(...)   REQUIRE(    (__VA_ARGS__))
#define    TRY_(...)   SECTION(    (__VA_ARGS__))
#define    EST_(...) BENCHMARK(    (__VA_ARGS__))

#define TAG_(...) TEST_CASE(__FILE__ ":" XTAL_S1_(__LINE__), TAG_N_(__VA_ARGS__))
#define TAG_N_(...)                                 XTAL_F1_(TAG_1_,__VA_ARGS__)
#define TAG_1_(...)                                             "[" __VA_ARGS__ "]"


namespace xtal
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///\returns the result of applying `std::complex` to the given arguments, \
zipping them together if vectorized. \

template <class X>
XTAL_FN2 complexion_f(X &&x)
XTAL_0EX
{
	using W = _std::complex<X>;
	return W(XTAL_REF_(x));
}
template <class X, class Y>
XTAL_FN2 complexion_f(X &&x, Y &&y)
XTAL_0EX
{
	using W = _std::complex<common_t<X, Y>>;
	return W(XTAL_REF_(x), XTAL_REF_(y));
}

//#if __has_include(<Eigen/Core>)
template <eigenvalue_q X>
XTAL_FN2 complexion_f(X &&x)
XTAL_0EX
{
	using W = _std::complex<eigenvalue_t<X>>;
	return flex_f<[] XTAL_1FN_(W)>(XTAL_REF_(x));
}
template <eigenvalue_q X, eigenvalue_q Y>
XTAL_FN2 complexion_f(X &&x, Y &&y)
XTAL_0EX
{
	using W = _std::complex<eigenvalue_t<X, Y>>;
	return flex_f<[] XTAL_1FN_(W)>(XTAL_REF_(x), XTAL_REF_(y));
}
//#endif

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////




XTAL_ENV_(push)
namespace xtal::math::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

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
