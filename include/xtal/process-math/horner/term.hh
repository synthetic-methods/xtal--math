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
XTAL_LET term_f(W &&w, X &&x, Xs &&...xs)
XTAL_0EX
{
	using V = based_t<devolved_u<X>>;
	auto constexpr n_sign = V{N_sign};
	auto const     x_sign =  (n_sign *...* XTAL_REF_(xs));

	XTAL_IF0
	XTAL_0IF (simplex_number_q<W, X> and complex_number_q<decltype(x_sign)>) {
		auto const &[xs_re, xs_im] = involved_f(x_sign);
		auto const y_im = x*xs_im;
		auto const y_re = term_f(XTAL_REF_(w), XTAL_REF_(x), xs_re);
		return _std::complex{y_re, y_im};
	}
	XTAL_0IF (complex_number_q<W> and simplex_number_q<X, decltype(x_sign)>) {
		auto const y_im = w.imag();
		auto const y_re = term_f(XTAL_REF_(w).real(), XTAL_REF_(x), x_sign);
		return _std::complex{y_re, y_im};
	}
	XTAL_0IF_(dynamic) {
		using _std::fma;
		if constexpr (bond::operate<V>::use_FMA() and requires {fma(x_sign, x, w);}) {
			return fma(x_sign, XTAL_REF_(x), XTAL_REF_(w));
		}
		else {
			return x_sign*XTAL_REF_(x) + XTAL_REF_(w);
		}
	}
	XTAL_0IF_(default) {
		return x_sign*XTAL_REF_(x) + XTAL_REF_(w);
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
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline)
		XTAL_SET function( auto &&...oo)
		XTAL_0EX -> decltype(auto)
		{
			return term_f(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
