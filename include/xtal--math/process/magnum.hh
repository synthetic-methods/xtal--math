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
struct magnum
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(_std::unsigned_integral auto const &o)
		noexcept -> auto
		{
			return o;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(_std::  signed_integral auto const &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			auto const v = o >> _fit::sign.shift;
			return (o^v) - v;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(real_variable_q auto const &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			return _xtd::copysign(o, _fit::alpha_1);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(complex_variable_q auto const &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			return root_f<2>(dot_f(o));
		}

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(real_variable_q auto &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			auto const o_sgn = _xtd::copysign(_fit::alpha_1, o);
			auto const o_mgn = o*o_sgn;
			o =    o_sgn;
			return o_mgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(complex_variable_q auto &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			auto [u, v] = dots_f<2>(o); o *= v; return u;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
using   magnum_t = process::confined_t<magnum<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(return,inline,let)
magnum_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return magnum_t<Ms...>::method_f(XTAL_REF_(oo)...);
	return magnum_t<>::template method_f<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
