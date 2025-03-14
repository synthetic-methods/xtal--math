#pragma once
#include "./any.hh"

#include "./tangent.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines `Tan[Pi #] &` and `Tanh[Pi #] &`, and their inverses.

Provides a counterpart to `pade::tangy` that is accurate in the limit for `Tanh`.

\tparam  M_ism
Specifies the underlying morphism \f$\in {1, 2}\f$,
generating either the circular or hyperbolic tangent.
*/

/*!
\brief   Produces the base-two tangy and its inverse.
\todo    Either define `subdivision` to apply `1/N_ind` prescaling, or supply an `attach`ed parameter to do so...
*/
template <int M_ism=0, int M_car=0>
struct  tangy;


////////////////////////////////////////////////////////////////////////////////

template <int M_iso, int M_car> requires (0 < M_iso)
struct tangy<M_iso, M_car>
:	process::lift<void
	,	tangent<M_iso, M_car>
	,	dilate<[] XTAL_1FN_(to) (one/bond::fit<>::patio_1)>
	>
{
};
template <int M_iso, int M_car> requires (M_iso < 0)
struct tangy<M_iso, M_car>
:	process::lift<void
	,	dilate<[] XTAL_1FN_(to) (one*bond::fit<>::patio_1)>
	,	tangent<M_iso, M_car>
	>
{
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=0, int M_car=0>
using tangy_t = process::confined_t<tangy<M_ism, M_car>>;

template <int M_ism=0, int M_car=0, int ...Ns>
XTAL_DEF_(let)
tangy_f = [] XTAL_1FN_(call) (tangy_t<M_ism, M_car>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
