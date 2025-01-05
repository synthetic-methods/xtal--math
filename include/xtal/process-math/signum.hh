#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "./dots.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct signum
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			return (0 < o) - (o < 0) + (o == 0);
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_number_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			return _xtd::copysign(_op::alpha_1, o);
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_number_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			return o*dot_f<-2>(o);
		}

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET edit(real_number_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			auto const o_sgn = function<Ns...>(o); o *= o_sgn; return o_sgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET edit(complex_number_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			auto [u, v] = dots_f<2>(o);
			auto const o_sgn = o*v;
			auto const o_mgn = XTAL_ALL_(o){u};
			o =    o_mgn;
			return o_sgn;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
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
template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET signum_e(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return signum_t<Ms...>::function(XTAL_REF_(oo)...);
	return signum_t<>::template     edit<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
