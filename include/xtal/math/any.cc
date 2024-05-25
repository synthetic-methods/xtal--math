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

/*/

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
	return zop_f<W>(XTAL_REF_(x));
}
template <eigenvalue_q X, eigenvalue_q Y>
XTAL_FN2 complexion_f(X &&x, Y &&y)
XTAL_0EX
{
	using W = _std::complex<eigenvalue_t<X, Y>>;
	return zop_f<W>(XTAL_REF_(x), XTAL_REF_(y));
}
/*/


template <template <class> class Y, class ...Xs>
XTAL_DEF_(return,inline)
XTAL_FN1 construxion_f(Xs &&...xs)
XTAL_0EX
{
	using W = common_t<Xs...>;
	return zop_f<Y<W>>(XTAL_REF_(xs)...);
}
template <template <class> class Y, class ...Xs> requires eigenvalue_q<Xs...>
XTAL_DEF_(return,inline)
XTAL_FN1 construxion_f(Xs &&...xs)
XTAL_0EX
{
	using W = common_t<eigenvalue_t<Xs>...>;
	return zop_f<Y<W>>(XTAL_REF_(xs)...);
}

XTAL_DEF_(return,inline)
XTAL_FN1 complexion_f(auto &&...xs)
XTAL_0EX
{
	return construxion_f<_std::complex>(XTAL_REF_(xs)...);
}
/***/


//template <class ...Xs> requires     common_q<Xs...> XTAL_FN2 complexion_f(Xs &&...xs) XTAL_0EX {return zop_f<_std::complex<    common_t<Xs...>>>(XTAL_REF_(xs)...);}
//template <class ...Xs> requires eigenvalue_q<Xs...> XTAL_FN2 complexion_f(Xs &&...xs) XTAL_0EX {return zop_f<_std::complex<eigenvalue_t<Xs...>>>(XTAL_REF_(xs)...);}
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
