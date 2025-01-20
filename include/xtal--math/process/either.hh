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
XTAL_DEF_(short)
XTAL_LET either_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return either_t<M_side>::template static_method<Ns...>(XTAL_REF_(oo)...);
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
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto const &x)
		noexcept -> auto
		{
			return static_method<Ns...>(x, XTAL_ALL_(x){});
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto const &x, auto const &y)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (M_side <= 0) return _std::min(x, y);
			XTAL_0IF (M_side == 1) return _std::max(x, y);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
