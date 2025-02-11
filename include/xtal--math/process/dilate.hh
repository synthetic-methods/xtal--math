#pragma once
#include "./any.hh"

#include "./root.hh"
#include "./signum.hh"
#include "./magnum.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <auto M_val=1> struct   dilate;
template <auto M_val=1> using    dilate_t = process::confined_t<dilate<M_val>>;
template <auto M_val=1>
XTAL_DEF_(return,inline,let)
dilate_f(auto &&o)
noexcept -> decltype(auto)
{
	return dilate_t<M_val>::method_f(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <auto M_val>
struct dilate
{
	XTAL_DEF_(return,inline,set)
	after_f(auto &&o)
	noexcept -> decltype(auto)
	{
		using _fit = bond::fit<decltype(o)>;
		auto constexpr n_val =   _fit::alpha_f(M_val);
		auto constexpr u     =       magnum_f(n_val);
		auto constexpr v     = (int) signum_f(n_val);
		return XTAL_REF_(o)*root_f<-v>(u);
	};

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			return after_f(XTAL_REF_(o));
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
