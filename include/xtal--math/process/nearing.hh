#pragma once
#include "./any.hh"

#include "./near.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
struct nearing
{
	using superkind = near<Ms...>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto const &u)
		noexcept -> XTAL_ALL_(u)
		{
			return u - S_::template static_method<Ns...>(u);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
using    nearing_t = process::confined_t<nearing<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET nearing_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return nearing_t<Ms...>::static_method(XTAL_REF_(oo)...);
	return nearing_t<>::template static_method<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
