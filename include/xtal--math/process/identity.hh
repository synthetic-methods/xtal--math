#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Provides the generic type signature for all `process`es. \

template <auto ...Ms> struct   identity;
template <auto ...Ms> using    identity_t = process::confined_t<identity<Ms...>>;


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct identity
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto &&...oo)
		noexcept -> decltype(auto)
		{
			return S_::template static_method<Ns...>(XTAL_REF_(oo)...);
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
