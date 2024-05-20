#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::math::horner
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_sign=1> struct term;
template <int N_sign=1> using  term_t = process::confined_t<term<N_sign>>;
template <int N_sign=1, additive_group_q W, multiplicative_group_q X, multiplicative_group_q ...Xs>
XTAL_DEF_(return,inline)
XTAL_LET term_f(W &&w, X &&x, Xs &&...xs)
XTAL_0EX -> decltype(auto)
{
	using _std::fma;

	using Xn = bond::seek_back_t<Xs...>;
	using re = bond::realize<Xn>;

//	based_t<Xn> const _s = N_sign; 

	if constexpr (re::N_fused and requires {fma((xs *...* N_sign), (x), w);}) {
		if (not _std::is_constant_evaluated()) {
			return fma((XTAL_REF_(xs) *...* N_sign), XTAL_REF_(x), XTAL_REF_(w));
		}
		else {
			return (XTAL_REF_(xs) *...* (N_sign*XTAL_REF_(x))) + XTAL_REF_(w);
		}
	}
	else {
		return (XTAL_REF_(xs) *...* (N_sign*XTAL_REF_(x))) + XTAL_REF_(w);
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
		XTAL_FN2 function( auto &&...oo)
		XTAL_0EX
		{
			return term_f(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
