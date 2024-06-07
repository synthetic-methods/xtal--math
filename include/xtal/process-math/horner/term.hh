#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::horner
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_sign=1> XTAL_TYP term;
template <int N_sign=1> XTAL_USE term_t = process::confined_t<term<N_sign>>;
template <int N_sign=1, additive_group_q W, multiplicative_group_q X, multiplicative_group_q ...Xs>
XTAL_DEF_(return,inline)
XTAL_FN1 term_f(W &&w, X &&x, Xs &&...xs)
XTAL_0EX
{
	using _std::fma;

	using Y = devolved_t<X>;
	using _op = bond::operate<Y>;
	based_t<Y> constexpr n_sign = N_sign;

	XTAL_IF0
	XTAL_0IF (none_q<Xs...>) {
		return XTAL_REF_(w) + n_sign*XTAL_REF_(x);
	}
	XTAL_0IF_(dynamic) {
		if constexpr (_op::use_FMA() and requires {fma((xs *...* n_sign), x, w);}) {
			return fma((XTAL_REF_(xs) *...* n_sign), XTAL_REF_(x), XTAL_REF_(w));
		}
		else {
			return (XTAL_REF_(xs) *...* (n_sign*XTAL_REF_(x))) + XTAL_REF_(w);
		}
	}
	XTAL_0IF_(default) {
		return (XTAL_REF_(xs) *...* (n_sign*XTAL_REF_(x))) + XTAL_REF_(w);
	}
};


////////////////////////////////////////////////////////////////////////////////
///\
Evaluates the term `(w + (x*xs...))` (using fused multiply-add, if supported by the compiler). \

///\
Used to define geometric/exponential series recursively via `(a[0] + b[0]*x*(...))`, \
starting from the kernel `a[N_limit]`. \

///\note\
Co/domain scaling can be effected by multiplying `a`/`b`, respectively. \

template <int N_sign>
struct term
{
//	static_assert(xtal::sign_p<N_sign, 1>);

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function( auto &&...oo)
		XTAL_0EX
		{
			return term_f(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
