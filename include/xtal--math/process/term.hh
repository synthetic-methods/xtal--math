#pragma once
#include "./any.hh"
#include "../atom/dot.hh"
#include "./square.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int M_alt=1, int M_pow=1, class W>
XTAL_DEF_(return,inline,let)
term_f(W &&w)
noexcept -> auto
{
	return XTAL_REF_(w);
}
template <int M_alt=1, int M_pow=1, class W, class X, class ...Xs>
XTAL_DEF_(return,inline,let)
term_f(W &&w, X &&x, Xs &&...xs)
noexcept -> auto
{
	using _xtd::plus_multiplies_f;
	auto constexpr then = [] XTAL_1FN_(call) (term_f<M_alt, M_pow>);
	using Y = unstruct_t<X, Xs...>;// NOTE: Constants interpreted as scalar quantities...
	XTAL_IF0

//	Decay integral/constants...
	XTAL_0IF (instant_q<W, X> and real_variable_q<Y>) {
		return then(static_cast<Y>(XTAL_REF_(w)), static_cast<Y>(XTAL_REF_(x)), XTAL_REF_(xs)...);
	}
	XTAL_0IF (instant_q<W   > and real_variable_q<Y>) {
		return then(static_cast<Y>(XTAL_REF_(w)),               (XTAL_REF_(x)), XTAL_REF_(xs)...);
	}
	XTAL_0IF (instant_q<   X> and real_variable_q<Y>) {
		return then(              (XTAL_REF_(w)), static_cast<Y>(XTAL_REF_(x)), XTAL_REF_(xs)...);
	}

//	Map groupoids...
	XTAL_0IF (atom::groupoid_q<W> and not atom::math::dot_q<W>) {
		return based_t<W>::template zip_from<then>(XTAL_REF_(w), XTAL_REF_(x), XTAL_REF_(xs)...);
	}
	XTAL_0IF (atom::groupoid_q<X> and not atom::math::dot_q<X>) {
		return based_t<X>::template zip_from<then>(XTAL_REF_(w), XTAL_REF_(x), XTAL_REF_(xs)...);
	}

//	Resolve parameters...
	XTAL_0IF (M_alt == 0) {
		return XTAL_REF_(w);
	}

	XTAL_0IF (M_alt ==  1 and M_pow == 1 and 0 == sizeof...(xs)) {
		return XTAL_REF_(w) + XTAL_REF_(x);
	}
	XTAL_0IF (M_alt == -1 and M_pow == 1 and 0 == sizeof...(xs)) {
		return XTAL_REF_(w) - XTAL_REF_(x);
	}
	XTAL_0IF (M_alt ==  1 and M_pow == 2 and 0 == sizeof...(xs)) {
		return  plus_multiplies_f( XTAL_REF_(w), x, x);
	}
	XTAL_0IF (M_alt == -1 and M_pow == 2 and 0 == sizeof...(xs)) {
		return -plus_multiplies_f(-XTAL_REF_(w), x, x);
	}
	/**/
	XTAL_0IF (M_pow == 2) {
		return then(XTAL_REF_(w), (XTAL_REF_(x) *...* XTAL_REF_(xs)));
	}
	/*/
	XTAL_0IF (M_alt ==  1 and M_pow == 2) {
		return XTAL_REF_(w) + square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));
	}
	XTAL_0IF (M_alt == -1 and M_pow == 2) {
		return XTAL_REF_(w) - square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));
	}
	/***/

	XTAL_0IF (M_alt ==  1 and M_pow == 1) {
		return  plus_multiplies_f( XTAL_REF_(w), XTAL_REF_(x), XTAL_REF_(xs)...);
	}
	XTAL_0IF (M_alt == -1 and M_pow == 1) {
		return -plus_multiplies_f(-XTAL_REF_(w), XTAL_REF_(x), XTAL_REF_(xs)...);
	}

	XTAL_0IF_(void)
};


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Evaluates the term `(w + (x*xs...))` (using fused multiply-add, if supported by the compiler).

Used to define geometric/exponential series recursively via `(a[0] + b[0]*x*(...))`,
starting from the kernel `a[M_limit]`.

\note    Co/domain scaling can be effected by multiplying `a`/`b`, respectively.
*/
template <auto ...Ms>
struct  term
{
//	static_assert(in_v<M_alt, 1,-1>);

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
using   term_t = process::confined_t<term<Ms...>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
