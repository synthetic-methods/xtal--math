#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
XTAL_TYP_(new) rate;

template <auto ...Ms>
XTAL_TYP_(let) rate_t = process::confined_t<rate<Ms...>>;

template <auto ...Ms>
XTAL_VAL_(let) rate_f = [] XTAL_1FN_(call) (rate_t<Ms...>::method);


////////////////////////////////////////////////////////////////////////////////

template <int M_ism>
struct rate<M_ism>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE
		template <auto ...Ns>
		XTAL_VAL_(return,inline,set)
		method(auto &&o)
		noexcept -> decltype(auto)
		{
			using U = XTAL_ALL_(o);
			XTAL_IF0
			XTAL_0IF (              real_q<U>) {return root_f<M_ism>(XTAL_REF_(o)        );}
			XTAL_0IF (           complex_q<U>) {return method<Ns...>(XTAL_REF_(o).real( ));}
			XTAL_0IF (atom::math::phason_q<U>) {return method<Ns...>(XTAL_REF_(o)     (1));}
			XTAL_0IF_(void)
		}
		template <auto ...Ns>
		XTAL_VAL_(return,inline,set)
		method(auto &&...oo)
		noexcept -> decltype(auto)
		{
			return root_f<M_ism>((one *...* rate_f<1>(XTAL_REF_(oo))));
		}

	};
};
template <>
struct rate<> : rate<1>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
