#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_alt=1, int M_pow=1, class W, multiplicative_group_q X, multiplicative_group_q ...Xs>
XTAL_DEF_(short)
XTAL_LET term_f(W &&w, X &&x, Xs &&...xs)
noexcept -> auto
{
	using X_   = absolve_u<X, Xs...>;
	using X_op = bond::operate<X>;
	X_ constexpr x_{M_alt};

	if constexpr (thunk_p<W, X_> or integer_q<W> and real_number_q<X_>) {
		return term_f<M_alt, M_pow>(static_cast<X_>(XTAL_REF_(w)), XTAL_REF_(x), XTAL_REF_(xs)...);
	}
	else {
		XTAL_IF0
		XTAL_0IF (M_pow == 0) {
			return XTAL_REF_(w);
		}
		XTAL_0IF (M_pow == 1) {
			auto const y = (XTAL_REF_(xs) *...* (x_));
			return _xtd::accumulator(XTAL_REF_(w), XTAL_REF_(x), y);
		}
		XTAL_0IF (M_pow == 2) {
			auto const y = (XTAL_REF_(xs) *...* (x_*XTAL_REF_(x)));
			return _xtd::accumulator(XTAL_REF_(w), y, y);
		}
		XTAL_0IF_(void)
	}
};


////////////////////////////////////////////////////////////////////////////////
///\
Evaluates the term `(w + (x*xs...))` (using fused multiply-add, if supported by the compiler). \

///\
Used to define geometric/exponential series recursively via `(a[0] + b[0]*x*(...))`, \
starting from the kernel `a[M_limit]`. \

///\note\
Co/domain scaling can be effected by multiplying `a`/`b`, respectively. \

template <auto ...Ms>
struct   term
{
//	static_assert(in_q<M_alt, 1,-1>);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(short,static)
		XTAL_LET function( auto &&...oo)
		noexcept -> decltype(auto)
		{
			return term_f<Ms...>(XTAL_REF_(oo)...);
		}

	};
};
template <auto ...Ms>
using    term_t = process::confined_t<term<Ms...>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
