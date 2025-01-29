#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <auto N>
XTAL_DEF_(return,inline,let)
factorial_f()
noexcept -> auto
{
	auto constexpr N_depth = sizeof(N) << 3;
	XTAL_IF0
	XTAL_0IF (N_depth == 0x80) {static_assert(N <= 34);}
	XTAL_0IF (N_depth == 0x40) {static_assert(N <= 20);}
	XTAL_0IF (N_depth == 0x20) {static_assert(N <= 12);}
	XTAL_IF0
	XTAL_0IF (N <= 1) {return one;}
	XTAL_0IF (2 <= N) {return N*factorial_f<N - 1>();}
}
template <integral_constant_q U>
XTAL_DEF_(return,inline,let)
factorial_f(U o)
noexcept -> auto
{
	return constant_t<factorial_f<U::value>()>{};
}
template <integral_variable_q U>
XTAL_DEF_(return,inline,let)
factorial_f(U o)
noexcept -> auto
{
	auto n{o}; while (--o) n *= o; return n;
}


////////////////////////////////////////////////////////////////////////////////
///\
Evaluates the factorial `(w + (x*xs...))` (using fused multiply-add, if supported by the compiler). \

///\
Used to define geometric/exponential series recursively via `(a[0] + b[0]*x*(...))`, \
starting from the kernel `a[M_limit]`. \

///\note\
Co/domain scaling can be effected by multiplying `a`/`b`, respectively. \

template <auto ...Ms>
struct   factorial
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		static_method(integral_variable_q auto &&o)
		noexcept -> decltype(auto)
		{
			return factorial_f(XTAL_REF_(o));
		}
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		static_method(real_variable_q auto &&o)
		noexcept -> decltype(auto)
		{
			using U_fix = bond::fixture<decltype(o)>;
			using U_sigma = typename U_fix::sigma_type;
			using U_alpha = typename U_fix::alpha_type;
			return static_cast<U_alpha>(factorial_f(static_cast<U_sigma>(XTAL_REF_(o))));
		}

	};
};
template <auto ...Ms>
using    factorial_t = process::confined_t<factorial<Ms...>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
