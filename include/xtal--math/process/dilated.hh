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
	XTAL_SET N_val = constant_t<M_val>{};

	template <auto f>
	XTAL_DEF_(short,static)
	XTAL_LET around_f(auto &&o)
	noexcept -> decltype(auto)
	{
		using _op = bond::operate<decltype(o)>;
		XTAL_LET n_val =   _op::alpha_f(N_val);
		XTAL_LET u     =       magnum_f(n_val);
		XTAL_LET v     = (int) signum_f(n_val);
		return f(XTAL_REF_(o)*root_f<-v>(u))*root_f<+v>(u);
	};

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		requires         in_n<XTAL_TRY_(S ::template function<Is...>  (XTAL_REF_(o)))>
		{
			return around_f<XTAL_FUN_(S_::template function<Is...>)>(XTAL_REF_(o));
		};

		template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET method(auto &&o) const
		noexcept -> decltype(auto)
		requires un_n<XTAL_TRY_(S::         template function<Is...>  (XTAL_REF_(o)))>
		and XTAL_TRY_(XTAL_ANY_(S_ const &).template   method<Is...>  (XTAL_REF_(o)))
		{
			return around_f<XTAL_FUN_(S_::template   method<Is...>)>(XTAL_REF_(o));
		};
		template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET method(auto &&o)
		noexcept -> decltype(auto)
		requires un_n<XTAL_TRY_(S::         template function<Is...>  (XTAL_REF_(o)))>
		and XTAL_TRY_(XTAL_ANY_(S_       &).template   method<Is...>  (XTAL_REF_(o)))
		{
			return around_f<XTAL_FUN_(S_::template   method<Is...>)>(XTAL_REF_(o));
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
