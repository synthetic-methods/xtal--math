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
	return dilate_t<M_val>::static_method(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <auto M_val>
struct dilate
{
	static auto constexpr N_val = constant_n<M_val>;

	XTAL_DEF_(return,inline,set)
	after_f(auto &&o)
	noexcept -> decltype(auto)
	{
		using _fix = bond::fixture<decltype(o)>;
		auto constexpr n_val =   _fix::alpha_f(N_val);
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
		static_method(auto &&o)
		noexcept -> auto
		{
			return after_f(XTAL_REF_(o));
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
