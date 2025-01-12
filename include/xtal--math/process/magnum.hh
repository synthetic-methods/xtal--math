#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "./dots.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
struct magnum
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(_std::unsigned_integral auto const &o)
		noexcept -> auto
		{
			return o;
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(_std::  signed_integral auto const &o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			auto const v = o >> _op::sign.shift;
			return (o^v) - v;
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_variable_q auto const &o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			return _xtd::copysign(o, _op::alpha_1);
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_variable_q auto const &o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			return dot_f<2>(o);
		}

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET edit(real_variable_q auto &o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			auto const o_sgn = _xtd::copysign(_op::alpha_1, o);
			auto const o_mgn = o*o_sgn;
			o =    o_sgn;
			return o_mgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET edit(complex_variable_q auto &o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			auto [u, v] = dots_f<2>(o); o *= v; return u;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
using    magnum_t = process::confined_t<magnum<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET magnum_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return magnum_t<Ms...>::function(XTAL_REF_(oo)...);
	return magnum_t<>::template function<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
