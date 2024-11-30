#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_alt=1> XTAL_TYP term;
template <int N_alt=1> XTAL_USE term_t = process::confined_t<term<N_alt>>;

template <int N_alt=1, additive_group_q W, multiplicative_group_q X, multiplicative_group_q ...Xs>
XTAL_DEF_(return,inline)
XTAL_LET term_f(W &&w, X &&x, Xs &&...xs)
noexcept -> auto
{
	using V = devolved_u<W, X, Xs...>;
	using _op = bond::operate<V>;

	auto constexpr v =  V{N_alt};
	auto const     y = (v *...* XTAL_REF_(xs));
	using Y = XTAL_ALL_(y);

	XTAL_IF0
	XTAL_0IF_(static) {
		return y*XTAL_REF_(x) + XTAL_REF_(w);
	}
	XTAL_0IF (_op::use_FMA() and requires {
		_std::fma(y, x, w);
		requires is_q<XTAL_ALL_(_std::fma(y, x, w)), XTAL_ALL_(y*x + w)>;
	}) {
		return _std::fma(y, XTAL_REF_(x), XTAL_REF_(w));
	}

	XTAL_0IF (complex_number_q<W> and simplex_number_q<X, Y>) {
		auto const &[w_re, w_im] = involved_f(XTAL_REF_(w));
		auto const & z_im = w_im;
		auto const   z_re = term_f(w_re, XTAL_REF_(x), y);
		return _std::complex{z_re, z_im};
	}

	/**/
	XTAL_0IF (simplex_number_q<W, X> and complex_number_q<Y>) {
		auto const &[y_re, y_im] = involved_f(y);
		auto const   z_im = x*y_im;
		auto const   z_re = term_f(XTAL_REF_(w), XTAL_REF_(x), y_re);
		return _std::complex{z_re, z_im};
	}
	XTAL_0IF (simplex_number_q<W, Y> and complex_number_q<X>) {
		return term_f(XTAL_REF_(w), y, XTAL_REF_(x));
	}
	/***/
	/**/
	XTAL_0IF (complex_number_q<W, X> and simplex_number_q<Y>) {
		auto const &[w_re, w_im] = involved_f(XTAL_REF_(w));
		auto const &[x_re, x_im] = involved_f(XTAL_REF_(x));
		auto const   z_re = term_f(w_re, x_re, y);
		auto const   z_im = term_f(w_im, x_im, y);
		return _std::complex{z_re, z_im};
	}
	XTAL_0IF (complex_number_q<W, Y> and simplex_number_q<X>) {
		return term_f(XTAL_REF_(w), y, XTAL_REF_(x));
	}
	/***/

	/*/
	XTAL_0IF (complex_number_q<W, X, Y>) {
		auto const &[w_re, w_im] = involved_f(XTAL_REF_(w));
		auto const &[x_re, x_im] = involved_f(XTAL_REF_(x));
		auto const &[y_re, y_im] = involved_f(y);
		auto const   z_re = term_f(term_f(w_re, x_re, y_re),-x_im, y_im);
		auto const   z_im = term_f(term_f(w_im, x_im, y_re), x_re, y_im);
		return _std::complex{z_re, z_im};
	}
	/***/

	XTAL_0IF_(else) {
	// <5% invocations in test-suite, most being `real + complex*complex`...
		return y*XTAL_REF_(x) + XTAL_REF_(w);
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

template <int N_alt>
struct term
{
//	static_assert(xtal::sign_p<N_alt, 1>);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function( auto &&...oo)
		noexcept -> decltype(auto)
		{
			return term_f(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
