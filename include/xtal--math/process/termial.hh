#pragma once
#include "./any.hh"

#include "./square.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Evaluates polynomial in `x` with coefficients `ks...`

Defines geometric/exponential series recursively via `(a[0] + b[0]*x*(...))`,
starting from the kernel `a[N_limit]`. \

\note    Uses fused multiply-add, if supported by the compiler).
\note    Co/domain scaling can be effected by multiplying `a`/`b`, respectively.
\note    Not to be confused with the Dijkstra's function `Binomial[n + 1, 2]`.
*/
template <int M_sgn=1, int M_pow=1>	struct  termial;
template <int M_sgn=1, int M_pow=1>	using   termial_t = process::confined_t<termial<M_sgn, M_pow>>;
template <int M_sgn=1, int M_pow=1>
XTAL_DEF_(return,inline,let)
termial_f(auto &&x, auto const &k, auto const &...ks)
noexcept -> auto
{
	XTAL_IF0
	XTAL_0IF (0 == sizeof...(ks)) {
		return XTAL_REF_(k);
	}
	XTAL_0IF (1 == M_pow) {
		return term_f<M_sgn>(k, x, termial_f<M_sgn>(x, ks...));
	}
	XTAL_0IF (2 == M_pow) {
		return termial_f<M_sgn>(square_f(XTAL_REF_(x)), k, ks...);
	}
	XTAL_0IF_(terminate)
};


////////////////////////////////////////////////////////////////////////////////
template <int M_sgn, int M_pow>
struct termial
{
//	static_assert(in_n<M_sgn, 1,-1>);

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
			return termial_f<M_sgn, M_pow>(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
