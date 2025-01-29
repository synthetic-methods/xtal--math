#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Wraps the super-`function`, unscaling/scaling the domain/codomain by `2^$1 Pi*$2`. \

///\note\
Because it invokes the super-`function` directly, \
it must be applied via `{compose,confined}` (etc) rather than `process::{lift,link}`.

template <auto M_val>
struct dilated;


////////////////////////////////////////////////////////////////////////////////

template <auto M_val>
struct dilated
{
	static auto constexpr N_val = constant_t<M_val>{};

	template <auto f>
	XTAL_DEF_(return,inline,set)
	around_f(auto &&o)
	noexcept -> decltype(auto)
	{
		using _fix = bond::fixture<decltype(o)>;
		auto constexpr n_val =   _fix::alpha_f(N_val);
		auto constexpr u     =       magnum_f(n_val);
		auto constexpr v     = (int) signum_f(n_val);
		return f(XTAL_REF_(o)*root_f<-v>(u))*root_f<+v>(u);
	};

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		static_method(auto &&o)
		noexcept -> decltype(auto)
		requires      in_n<requires {S ::template static_method<Is...>  (XTAL_REF_(o));}>
		{
			return around_f<[] XTAL_0FN_(alias) (S_::template static_method<Is...>)>(XTAL_REF_(o));
		};

		template <auto ...Is>
		XTAL_DEF_(return,inline,let)
		method(auto &&o) const
		noexcept -> decltype(auto)
		requires      un_n<requires {S ::template static_method<Is...>  (XTAL_REF_(o));}>
		and requires (S_ const &s_) {s_ .template        method<Is...>  (XTAL_REF_(o));}
		{
			return around_f<[] XTAL_0FN_(alias) (S_::template        method<Is...>)>(XTAL_REF_(o));
		};
		template <auto ...Is>
		XTAL_DEF_(return,inline,let)
		method(auto &&o)
		noexcept -> decltype(auto)
		requires      un_n<requires {S ::template static_method<Is...>  (XTAL_REF_(o));}>
		and requires (S_       &s_) {s_ .template        method<Is...>  (XTAL_REF_(o));}
		{
			return around_f<[] XTAL_0FN_(alias) (S_::template        method<Is...>)>(XTAL_REF_(o));
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
