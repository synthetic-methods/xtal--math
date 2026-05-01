#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Evaluates the exponential `Times[oo]^x`.
*/
template <auto ...Ms>
struct  exponential;


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct  exponential
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method(auto x)
		noexcept -> auto
		{
			return one;
		}
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method(cardinal_variable_q auto x, auto &&...oo)
		noexcept -> auto
		{
			auto u = (one *...* XTAL_REF_(oo));
			XTAL_ALL_(u) w{one};
			while (x) {
				if (x&1U) {
					w *= u;
				}
				u  *= u;
				x >>= 1;
			}
			return w;
		}

	};
};
////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
XTAL_TYP_(let) exponential_t = process::confined_t<exponential<Ms...>>;

template <auto ...Ms>
XTAL_DEF_(let) exponential_f = [] XTAL_1FN_(call) (exponential_t<Ms...>::method);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
