#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_alt=1, int M_pow=1, class W, class X, class ...Xs>
XTAL_DEF_(return,inline,let)
term_f(W &&w, X &&x, Xs &&...xs)
noexcept -> auto
{
	if constexpr (0 == sizeof...(Xs)) {
		XTAL_IF0
		XTAL_0IF (M_pow == 0) {
			return XTAL_REF_(w);
		}
		XTAL_0IF (M_pow == 1) {
			return XTAL_REF_(w) + XTAL_REF_(x);
		}
		XTAL_0IF (M_pow == 2) {
			return term_f<M_alt>(XTAL_REF_(w), x, x);
		}
	}
	else {
		using Y = absolve_u<Xs...>;// NOTE: Constants interpreted as scalar quantities...
		XTAL_IF0
		XTAL_0IF (constant_q<W, X> or integral_variable_q<W, X> and real_variable_q<Y>) {
			return term_f<M_alt, M_pow>(static_cast<Y>(XTAL_REF_(w)), static_cast<Y>(XTAL_REF_(x)), XTAL_REF_(xs)...);
		}
		XTAL_0IF (constant_q<W   > or integral_variable_q<W   > and real_variable_q<Y>) {
			return term_f<M_alt, M_pow>(static_cast<Y>(XTAL_REF_(w)), XTAL_REF_(x), XTAL_REF_(xs)...);
		}
		XTAL_0IF (constant_q<   X> or integral_variable_q<   X> and real_variable_q<Y>) {
			return term_f<M_alt, M_pow>(XTAL_REF_(w), static_cast<Y>(XTAL_REF_(x)), XTAL_REF_(xs)...);
		}
		XTAL_0IF (M_pow == 0  or M_alt ==  0) {
			return XTAL_REF_(w);
		}
		XTAL_0IF (M_pow == 1 and M_alt ==  1) {
			return accumulator_f(XTAL_REF_(w), XTAL_REF_(x), XTAL_REF_(xs)...);
		}
		XTAL_0IF (M_pow == 1 and M_alt == -1) {
			return accumulator_f(XTAL_REF_(w), XTAL_REF_(x), XTAL_REF_(xs)..., Y{M_alt});
		}
		XTAL_0IF (M_pow == 2 and M_alt ==  1) {
			auto const y = (XTAL_REF_(xs) *...* XTAL_REF_(x));
			return accumulator_f(XTAL_REF_(w), y, y);
		}
		XTAL_0IF (M_pow == 2 and M_alt == -1) {
			auto const y = (XTAL_REF_(xs) *...* XTAL_REF_(x));
			return accumulator_f(XTAL_REF_(w), y, y, Y{M_alt});
		}
		XTAL_0IF_(terminate)
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
//	static_assert(in_n<M_alt, 1,-1>);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&...oo)
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
