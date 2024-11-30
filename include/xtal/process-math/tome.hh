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


template <int N_sign=1> XTAL_TYP tome;
template <int N_sign=1> XTAL_USE tome_t = process::confined_t<tome<N_sign>>;
template <int N_sign=1>
XTAL_DEF_(return,inline)
XTAL_LET tome_f(auto &&w, auto &&k, auto &&...ks)
noexcept -> decltype(auto)
{
	if constexpr (0 == sizeof...(ks)) {
		return static_cast<XTAL_ALL_(k)>(XTAL_REF_(k));
	}
	else {
		return term_f<N_sign>(k, w, tome_f<N_sign>(w, ks...));
	}
};


////////////////////////////////////////////////////////////////////////////////
template <int N_sign>
struct tome
{
//	static_assert(xtal::sign_p<N_sign, 1>);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&...oo)
		noexcept -> decltype(auto)
		{
			return tome_f<N_sign>(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
