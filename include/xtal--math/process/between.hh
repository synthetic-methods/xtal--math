#pragma once
#include "./any.hh"

#include "./magnum.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///\returns the argument restricted to the closed interval `[M_dn, M_up]`. \

template <int M_dn=0, int M_up=0>
XTAL_DEF_(return,inline,let)
between_f(auto &&u)
noexcept -> auto
{
	XTAL_IF0
	XTAL_0IF_(consteval) {
		return u < M_dn? M_dn: M_up < u? M_up: u;
	}
	XTAL_0IF_(else) {
		using U = XTAL_ALL_(u);
		U constexpr  n0 = M_dn;
		U constexpr  n1 = M_up;
		U constexpr n01 = M_dn + M_up;
		if constexpr (integral_variable_q<U>) {
			return (magnum_f(u - n0) - magnum_f(u - n1) + n01) >> 1;
		}
		else {
			return (magnum_f(u - n0) - magnum_f(u - n1) + n01)*0.5f;
		}
	}
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct   between
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
			return between_f<Ms...>(XTAL_REF_(oo)...);
		}

	};
};
template <auto ...Ms>
using    between_t = process::confined_t<between<Ms...>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
