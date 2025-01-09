#pragma once
#include "./any.hh"

#include "./root.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Combines `power` (when `0 < M_ism`) and `root` (when `M_ism < 0`) as mutual inverses, \
equivalent to `#^m^Sgn@m &`. \


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, auto ...Ms>
struct   nomial;

template <int M_ism, auto ...Ms> requires (M_ism <  0) struct nomial<M_ism, Ms...> : root<-M_ism, Ms...> {};
template <int M_ism, auto ...Ms> requires (0 <= M_ism) struct nomial<M_ism, Ms...> : power<M_ism, Ms...> {};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
using    nomial_t = process::confined_t<nomial<Ms...>>;

template <int ...Ms>
XTAL_DEF_(short)
XTAL_LET nomial_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return nomial_t<Ms...>::function(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
