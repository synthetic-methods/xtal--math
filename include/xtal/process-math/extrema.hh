#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_side=0>
struct   extrema;

template <int M_side=0>
using    extrema_t = process::confined_t<extrema<M_side>>;

template <int M_side=0, auto ...Ns>
XTAL_DEF_(short)
XTAL_LET extrema_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return extrema_t<M_side>::template function<Ns...>(XTAL_REF_(oo)...);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_side>
struct extrema
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &x)
		noexcept -> auto
		{
			return function<Ns...>(x, XTAL_ALL_(x){});
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &x, auto const &y)
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
