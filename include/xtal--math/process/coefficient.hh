#pragma once
#include "./any.hh"

#include "./term.hh"
#include "./square.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Multiplies the leading argument with the result of the parent-`method` applied to the trailing arguments.
\todo    Implement positional parameter `M_pos`?
*/
template <auto ...Ms>
struct coefficient;


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct coefficient
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&o, auto &&...oo)
		noexcept -> auto
		{
			return S_::template method<Ns...>(XTAL_REF_(oo)...)*XTAL_REF_(o);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
