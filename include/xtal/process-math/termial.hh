#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Evaluates polynomial in `w` with coefficients `ks...` \
(using fused multiply-add, if supported by the compiler). \

///\
Used to define geometric/exponential series recursively via `(a[0] + b[0]*x*(...))`, \
starting from the kernel `a[N_limit]`. \

///\note\
Co/domain scaling can be effected by multiplying `a`/`b`, respectively. \

///\note\
Not to be confused with the Dijkstra's function `Binomial[n + 1, 2]`. \


template <int M_sgn=1>	struct   termial;
template <int M_sgn=1>	using    termial_t = process::confined_t<termial<M_sgn>>;
template <int M_sgn=1>
XTAL_DEF_(short)
XTAL_LET termial_f(auto &&w, auto &&k, auto &&...ks)
noexcept -> decltype(auto)
{
	if constexpr (0 == sizeof...(ks)) {
		return static_cast<XTAL_ALL_(k)>(XTAL_REF_(k));
	}
	else {
		return term_f<M_sgn>(k, w, termial_f<M_sgn>(w, ks...));
	}
};


////////////////////////////////////////////////////////////////////////////////
template <int M_sgn>
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
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&...oo)
		noexcept -> decltype(auto)
		{
			return termial_f<M_sgn>(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
