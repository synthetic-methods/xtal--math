#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Evaluates the monomial `Module[{ms = ##}, Times@@({##}^{ms}) &] &`.
*/
template <int ...Ms> XTAL_TYP_(new)  monomial;
template <int ...Ms> XTAL_TYP_(let)  monomial_t = process::confined_t<monomial<Ms...>>;
template <int ...Ms> requires (0 == sizeof...(Ms))
XTAL_DEF_(return,inline,set)
monomial_f()
noexcept -> auto
{
	return one;
}
template <int M> requires (0 == M)
XTAL_DEF_(return,inline,set)
monomial_f(objective_q auto &&o)
noexcept -> XTAL_ALL_(o)
{
	return XTAL_ALL_(o) {one};
}
template <int M> requires (1 == M)
XTAL_DEF_(return,inline,set)
monomial_f(objective_q auto &&o)
noexcept -> XTAL_ALL_(o)
{
	return XTAL_REF_(o);
}
template <int M> requires (2 == M)
XTAL_DEF_(return,inline,set)
monomial_f(objective_q auto &&o)
noexcept -> XTAL_ALL_(o)
{
	using _fit = bond::fit<decltype(o)>;
	if constexpr (complex_variable_q<decltype(o)>) {
		auto const x  = _std::real(o), xx = monomial_f<2>(x);
		auto const y  = _std::imag(o), yy = monomial_f<2>(y);
		auto const u  = xx - yy;
		auto const v  = two*x*y;
		return {u, v};
	}
	else {
		return o*o;
	}
}
template <int M> requires (3 == M)
XTAL_DEF_(return,inline,set)
monomial_f(objective_q auto &&o)
noexcept -> XTAL_ALL_(o)
{
	using _fit = bond::fit<decltype(o)>;
	if constexpr (complex_variable_q<decltype(o)>) {
		using U = XTAL_ALL_(o);
		using V = typename U::value_type;
		auto const x  = _std::real(o), xx = monomial_f<2>(x);
		auto const y  = _std::imag(o), yy = monomial_f<2>(y);
		auto const u  = _xtd::plus_multiplies_f(xx, yy, V(-3))* x;
		auto const v  = _xtd::plus_multiplies_f(yy, xx, V(-3))*-y;
		return {u, v};
	}
	else {
		return o*o*o;
	}
}
template <int M> requires (4 <= M)
XTAL_DEF_(return,inline,set)
monomial_f(objective_q auto &&o)
noexcept -> XTAL_ALL_(o)
{
//	TODO: Loop/unroll modulo `Floor[(32 - 1)/Log[{2,3}]]`?
	XTAL_IF0
	XTAL_0IF (0 == M%3) {return monomial_f<M / 3>(monomial_f<3>(XTAL_REF_(o)));}
	XTAL_0IF (0 == M%2) {return monomial_f<M / 2>(monomial_f<2>(XTAL_REF_(o)));}
	XTAL_0IF_(else)     {return monomial_f<M - 1>(o)*(o);}
}

template <int M> requires (0 <= M)
XTAL_DEF_(return,inline,set)
monomial_f(subjective_q auto &&o)
noexcept -> auto
{
	return monomial_f(objective_f(XTAL_REF_(o)));
}

template <int M> requires (M <  0)
XTAL_DEF_(return,inline,set)
monomial_f(auto &&o)
noexcept -> auto
{
	return monomial_f<-M>(one/XTAL_REF_(o));
}

template <int M, int ...Ms>
XTAL_DEF_(return,inline,set)
monomial_f(auto &&o, auto &&...oo)
noexcept -> auto
requires (1 <= sizeof...(Ms)) and (1 <= sizeof...(oo))
{
	auto constexpr L = sizeof...(Ms);
	static_assert(L == sizeof...(oo));
	XTAL_IF0
	XTAL_0IF (0 == L) {return monomial_f<M>(XTAL_REF_(o));}
	XTAL_0IF (1 <= L) {return monomial_f<M>(XTAL_REF_(o))*monomial_f<Ms...>(XTAL_REF_(oo)...);}
}



////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
struct  monomial
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	public:

		template <int ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&...oo)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (0 == sizeof...(Ms)) {return monomial_f<Ns...>(XTAL_REF_(oo)...);}
			XTAL_0IF (0 == sizeof...(Ns)) {return monomial_f<Ms...>(XTAL_REF_(oo)...);}
			XTAL_0IF_(else)               {return one;}
			
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
