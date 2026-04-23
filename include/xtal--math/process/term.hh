#pragma once
#include "./any.hh"
#include "../atom/dot.hh"
#include "./square.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int M_add=1, int M_mul=1, class W>
XTAL_DEF_(return,inline,let)
term_f(W &&w)
noexcept -> auto
{
	return XTAL_REF_(w);
}
template <int M_add=1, int M_mul=1, class W, class X, class ...Xs>
XTAL_DEF_(return,inline,let)
term_f(W &&w, X &&x, Xs &&...xs)
noexcept -> auto
{
	static_assert(-2 <= M_add and M_add <= +2);
	static_assert(-2 <= M_mul and M_mul <= +2);

	using _xtd::plus_multiplies_f;
	auto constexpr then = [] XTAL_1FN_(call) (term_f<M_add, M_mul>);
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
	XTAL_0IF (M_mul == 0 or same_q<X, decltype(zero)>) {
		XTAL_IF0
		XTAL_0IF (M_add ==  1) {return          (XTAL_REF_(w));}
		XTAL_0IF (M_add ==  2) {return  square_f(XTAL_REF_(w));}
		XTAL_0IF (M_add == -1) {return -        (XTAL_REF_(w));}
		XTAL_0IF (M_add == -2) {return -square_f(XTAL_REF_(w));}
	}
	XTAL_0IF (M_add == 0 or same_q<W, decltype(zero)>) {
		XTAL_IF0
		XTAL_0IF (M_mul ==  1) {return          ((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (M_mul ==  2) {return  square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (M_mul == -1) {return -        ((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (M_mul == -2) {return -square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
	}
	XTAL_0IF (0 == sizeof...(Xs)) {
		XTAL_IF0
		XTAL_0IF (M_add ==  1 and M_mul ==  1) {return   XTAL_REF_(w) + XTAL_REF_(x);}
		XTAL_0IF (M_add == -1 and M_mul == -1) {return - XTAL_REF_(x) - XTAL_REF_(w);}
		XTAL_0IF (M_add == -1 and M_mul ==  1) {return   XTAL_REF_(x) - XTAL_REF_(w);}
		XTAL_0IF (M_add ==  1 and M_mul == -1) {return   XTAL_REF_(w) - XTAL_REF_(x);}
		XTAL_0IF (M_add ==  1 and M_mul ==  2) {return   plus_multiplies_f(         (XTAL_REF_(w)), x, x);}
		XTAL_0IF (M_add == -1 and M_mul == -2) {return - plus_multiplies_f(         (XTAL_REF_(w)), x, x);}
		XTAL_0IF (M_add == -1 and M_mul ==  2) {return   plus_multiplies_f(-        (XTAL_REF_(w)), x, x);}
		XTAL_0IF (M_add ==  1 and M_mul == -2) {return - plus_multiplies_f(-        (XTAL_REF_(w)), x, x);}
		XTAL_0IF (M_add ==  2 and M_mul ==  2) {return   plus_multiplies_f( square_f(XTAL_REF_(w)), x, x);}
		XTAL_0IF (M_add == -2 and M_mul == -2) {return - plus_multiplies_f( square_f(XTAL_REF_(w)), x, x);}
		XTAL_0IF (M_add == -2 and M_mul ==  2) {return   plus_multiplies_f(-square_f(XTAL_REF_(w)), x, x);}
		XTAL_0IF (M_add ==  2 and M_mul == -2) {return - plus_multiplies_f(-square_f(XTAL_REF_(w)), x, x);}
	}
	XTAL_0IF (1 <= sizeof...(Xs)) {
		XTAL_IF0
		/**/
		XTAL_0IF (                M_mul ==  2) {return   then(XTAL_REF_(w),  (XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (                M_mul == -2) {return   then(XTAL_REF_(w), -(XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (M_add ==  2                ) {return   plus_multiplies_f( square_f(XTAL_REF_(w)), XTAL_REF_(xs)...);}
		XTAL_0IF (M_add == -2                ) {return   plus_multiplies_f(-square_f(XTAL_REF_(w)), XTAL_REF_(xs)...);}
		/*/
		XTAL_0IF (M_add ==  1 and M_mul ==  2) {return   XTAL_REF_(w) + square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (M_add == -1 and M_mul == -2) {return - XTAL_REF_(w) - square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (M_add == -1 and M_mul ==  2) {return - XTAL_REF_(w) + square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		XTAL_0IF (M_add ==  1 and M_mul == -2) {return   XTAL_REF_(w) - square_f((XTAL_REF_(x) *...* XTAL_REF_(xs)));}
		/***/

		XTAL_0IF (M_add ==  1 and M_mul ==  1) {return   plus_multiplies_f( XTAL_REF_(w),  XTAL_REF_(x), XTAL_REF_(xs)...);}
		XTAL_0IF (M_add == -1 and M_mul == -1) {return - plus_multiplies_f( XTAL_REF_(w),  XTAL_REF_(x), XTAL_REF_(xs)...);}
		XTAL_0IF (M_add == -1 and M_mul ==  1) {return   plus_multiplies_f(-XTAL_REF_(w),  XTAL_REF_(x), XTAL_REF_(xs)...);}
		XTAL_0IF (M_add ==  1 and M_mul == -1) {return   plus_multiplies_f( XTAL_REF_(w), -XTAL_REF_(x), XTAL_REF_(xs)...);}
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
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,let)
		method(auto &&...oo)
		const noexcept -> decltype(auto)
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
