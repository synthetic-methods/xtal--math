#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///\returns the argument modulo `M_lim`. \

template <int M_lim>
XTAL_DEF_(return,inline,let)
modulo_f(auto &&u)
noexcept -> auto
{
	return (((XTAL_REF_(u)%M_lim) + M_lim)%M_lim);
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct  modulo
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&...oo)
		noexcept -> decltype(auto)
		{
			return modulo_f<Ms...>(XTAL_REF_(oo)...);
		}

	};
};
template <auto ...Ms>
using   modulo_t = process::confined_t<modulo<Ms...>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
