#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*!
*/
template <                  class ..._s>	struct  dash;
template <         class S, class ..._s>	using   dash_s  = bond::compose_s<packed_t< _s...>, dash<S>>;
template <         class S, class ...Ts>	concept dash_p  = bond::tag_outer_p<dash, Ts...> and bond::tab_inner_p<S, Ts...>;
template <class T, class S             >	concept dash_q  = bond::tag_outer_p<dash, T    > and bond::tab_inner_p<S, T    >;
template <         class S             >	struct  dash<S> : bond::compose<bond::tag<dash>, bond::tab<S>> {};


////////////////////////////////////////////////////////////////////////////////

template <class S>
XTAL_DEF_(let) dash_f = [] (auto &&...oo) XTAL_0FN_(to) (dash_s<S, XTAL_ALL_(oo)...>(XTAL_REF_(oo)...));


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
