#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
struct signum
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u)
		noexcept -> XTAL_ALL_(u)
		{
			return (0 < u) - (u < 0) + (u == 0);
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_number_q auto const &u)
		noexcept -> XTAL_ALL_(u)
		{
			return _xtd::signbit(u);
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_number_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			return {function<Ns...>(o.real()), function<Ns...>(o.imag())};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
using    signum_t = process::confined_t<signum<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET signum_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return signum_t<Ms...>::function(XTAL_REF_(oo)...);
	return signum_t<>::template function<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
