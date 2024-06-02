#pragma once
#include "./any.hh"
#include "./term.hh"





XTAL_ENV_(push)
namespace xtal::process::math::horner
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_sign=1> struct polynomial;
template <int N_sign=1> using  polynomial_t = process::confined_t<polynomial<N_sign>>;
template <int N_sign=1>
XTAL_DEF_(return,inline)
XTAL_FN1 polynomial_f(auto &&w, auto &&k, auto &&...ks)
XTAL_0EX
{
	if constexpr (0 == sizeof...(ks)) {
		return k;
	}
	else {
		return term_f<N_sign>(k, w, polynomial_f<N_sign>(w, ks...));
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
struct polynomial
{
//	static_assert(xtal::sign_p<N_sign, 1>);

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_LET function = polynomial_f<N_sign>;

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
