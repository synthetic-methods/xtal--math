#pragma once
#include <xtal/any.hh>
#include <xtal/all.hh>
#include <Eigen/Dense>

#include "./etc.hh"


XTAL_ENV_(push)
namespace xtal::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

///\returns the result of applying the `f` to `...xs`, \
zipping them together if vectorized. \

///\note\
Provides experimental support for `Eigen` (via `\.(?:un|bin|tern)aryExpr`), \
but can be specialized to support additional/custom data-types. \

///\todo\
Restrict `eigenvalue_q` to `Array`-derived types.

template <class F, eigenvalue_q X0                                  > XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, X0 &&x0                  ) XTAL_0EX {return XTAL_REF_(x0).  unaryExpr(                              XTAL_REF_(f));}
template <class F, eigenvalue_q X0, eigenvalue_q X1                 > XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, X0 &&x0, X1 &&x1         ) XTAL_0EX {return XTAL_REF_(x0). binaryExpr(XTAL_REF_(x1),                XTAL_REF_(f));}
template <class F, eigenvalue_q X0, eigenvalue_q X1, eigenvalue_q X2> XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, X0 &&x0, X1 &&x1, X2 &&x2) XTAL_0EX {return XTAL_REF_(x0).ternaryExpr(XTAL_REF_(x1), XTAL_REF_(x2), XTAL_REF_(f));}
template <class F, class ...Xs> requires some_q<Xs...>                XTAL_DEF_(inline) XTAL_LET hoist_f(F &&f, Xs &&...xs)                XTAL_0EX {return             XTAL_REF_(f) (XTAL_REF_(xs)...);}
template <class F, class ...Xs> requires some_q<Xs...>                XTAL_DEF_(inline) XTAL_LET hoist_f(       Xs &&...xs)                XTAL_0EX {return      hoist_f(invoke_f<F>, XTAL_REF_(xs)...);}


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


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
