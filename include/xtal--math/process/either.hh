#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_side=0>
struct   either;

template <int M_side=0>
using    either_t = process::confined_t<either<M_side>>;

template <int M_side=0, auto ...Ns>
XTAL_DEF_(return,inline,let)
either_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return either_t<M_side>::template method_f<Ns...>(XTAL_REF_(oo)...);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_side>
struct either
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &x)
		noexcept -> auto
		{
			return method_f<Ns...>(x, XTAL_ALL_(x){});
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &x, auto const &y)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(x), decltype(y)>;
			XTAL_IF0
			XTAL_0IF (M_side <= 0) return _fit::minimum_f(x, y);
			XTAL_0IF (M_side == 1) return _fit::maximum_f(x, y);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
