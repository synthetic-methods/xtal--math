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

#if __has_include(<Eigen/Core>)
template <class         T >	XTAL_USE     eigenclass_t =	         Eigen::internal::traits<based_t<T>>;
template <class         T >	XTAL_USE     eigenvalue_t =	typename Eigen::internal::traits<based_t<T>>::Scalar;
template <class      ...Ts>	XTAL_ASK     eigenclass_q =	complete_q<eigenclass_t<Ts>...>;
template <class      ...Ts>	XTAL_ASK     eigenvalue_q =	complete_q<eigenvalue_t<Ts>...>;//TODO: Restrict to `Array`-derived.
#else
template <class         T >	XTAL_USE     eigenclass_t =	void;
template <class         T >	XTAL_USE     eigenvalue_t =	void;
template <class      ...Ts>	XTAL_ASK     eigenclass_q =	false;
template <class      ...Ts>	XTAL_ASK     eigenvalue_q =	false;
#endif

namespace _retail
{
	template <eigenvalue_q T>
	XTAL_TYP devolve<T> {using value_type = eigenvalue_t<T>;};

}

///\returns the result of applying the `f` to `...xs`, \
zipping them together if vectorized. \

///\note\
Provides experimental support for `Eigen` (via `\.(?:un|bin|tern)aryExpr`), \
but can be specialized to support additional/custom data-types. \

///\todo\
Restrict `eigenvalue_q` to `Array`-derived types.

#if __has_include(<Eigen/Core>)
template <class F, eigenvalue_q X0                                  > XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, X0 &&x0                  ) XTAL_0EX {return XTAL_REF_(x0).  unaryExpr(                              XTAL_REF_(f));}
template <class F, eigenvalue_q X0, eigenvalue_q X1                 > XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, X0 &&x0, X1 &&x1         ) XTAL_0EX {return XTAL_REF_(x0). binaryExpr(XTAL_REF_(x1),                XTAL_REF_(f));}
template <class F, eigenvalue_q X0, eigenvalue_q X1, eigenvalue_q X2> XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, X0 &&x0, X1 &&x1, X2 &&x2) XTAL_0EX {return XTAL_REF_(x0).ternaryExpr(XTAL_REF_(x1), XTAL_REF_(x2), XTAL_REF_(f));}
#endif
template <class F, class ...Xs> requires some_q<Xs...>                XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, Xs &&...xs)                XTAL_0EX {return                                       XTAL_REF_(f) (XTAL_REF_(xs)...);}
template <class F, class ...Xs> requires some_q<Xs...>                XTAL_DEF_(inline) XTAL_LET hoist_f(       Xs &&...xs)                XTAL_0EX {return                                  hoist_f(invoke_f<F>, XTAL_REF_(xs)...);}



template <template <class> class Y, class ...Xs>
XTAL_DEF_(return,inline)
XTAL_FN1 construxion_f(Xs &&...xs)
XTAL_0EX
{
	using W = common_t<Xs...>;
	return hoist_f<Y<W>>(XTAL_REF_(xs)...);
}
template <template <class> class Y, class ...Xs> requires eigenvalue_q<Xs...>
XTAL_DEF_(return,inline)
XTAL_FN1 construxion_f(Xs &&...xs)
XTAL_0EX
{
	using W = common_t<eigenvalue_t<Xs>...>;
	return hoist_f<Y<W>>(XTAL_REF_(xs)...);
}

XTAL_DEF_(return,inline)
XTAL_FN1 complexion_f(auto &&...xs)
XTAL_0EX
{
	return construxion_f<_std::complex>(XTAL_REF_(xs)...);
}


//template <class ...Xs> requires     common_q<Xs...> XTAL_FN2 complexion_f(Xs &&...xs) XTAL_0EX {return hoist_f<_std::complex<    common_t<Xs...>>>(XTAL_REF_(xs)...);}
//template <class ...Xs> requires eigenvalue_q<Xs...> XTAL_FN2 complexion_f(Xs &&...xs) XTAL_0EX {return hoist_f<_std::complex<eigenvalue_t<Xs...>>>(XTAL_REF_(xs)...);}
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
