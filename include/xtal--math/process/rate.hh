#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int M_ism=0, int M_cut=0>
struct rate
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(return,inline,set)
		static_method(simplex_field_q auto const &o_re, simplex_field_q auto const &o_im)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (M_ism == -1) {return o_im*root_f<-1, M_cut>(o_re);}
			XTAL_0IF (M_ism == +1) {return o_re*root_f<-1, M_cut>(o_im);}
			XTAL_0IF_(void)
		}
		template <int ...Ns>
		XTAL_DEF_(return,inline,set)
		static_method(complex_field_q auto const &o)
		noexcept -> auto
		{
			return static_method<Ns...>(o.real(), o.imag());
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
using    rate_t = process::confined_t<rate<Ms...>>;

template <auto ...Ms>
XTAL_DEF_(return,inline,let)
rate_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return rate_t<Ms...>::static_method(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
