#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_alt=1, additive_group_q W, multiplicative_group_q X, multiplicative_group_q ...Xs>
XTAL_DEF_(short)
XTAL_LET term_f(W &&w, X &&x, Xs &&...xs)
noexcept -> auto
{
	absolve_u<W, X, Xs...> constexpr v{N_alt};
	return _xtd::fam(XTAL_REF_(w), XTAL_REF_(x), (v *...* XTAL_REF_(xs)));
};


////////////////////////////////////////////////////////////////////////////////
///\
Evaluates the term `(w + (x*xs...))` (using fused multiply-add, if supported by the compiler). \

///\
Used to define geometric/exponential series recursively via `(a[0] + b[0]*x*(...))`, \
starting from the kernel `a[N_limit]`. \

///\note\
Co/domain scaling can be effected by multiplying `a`/`b`, respectively. \

template <int N_alt=1>
struct   term
{
//	static_assert(in_q<N_alt, 1,-1>);

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
			return term_f<N_alt>(XTAL_REF_(oo)...);
		}

	};
};
template <int N_alt=1>
using    term_t = process::confined_t<term<N_alt>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
