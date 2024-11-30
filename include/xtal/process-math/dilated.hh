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

template <int N_two=0, int N_two_pi=0>
struct dilated;


////////////////////////////////////////////////////////////////////////////////

template <int N_two, int N_two_pi>
struct dilated
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&o)
		noexcept -> decltype(auto)
			requires (not XTAL_TRY_(S::template function<Is...>(XTAL_REF_(o))))
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr n = _op::diplo_f(-N_two)*_op::template patio_f<-N_two_pi>(2, 1);
			auto constexpr u = _op::diplo_f(+N_two)*_op::template patio_f<+N_two_pi>(2, 1);
			
			return S_::template method<Is...>(XTAL_REF_(o)*n)*u;
		};
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&o) const
		noexcept -> decltype(auto)
			requires (not XTAL_TRY_(S::template function<Is...>(XTAL_REF_(o))))
			and XTAL_TRY_(XTAL_ANY_(S_ const &).template method<Is...>(XTAL_REF_(o)))
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr n = _op::diplo_f(-N_two)*_op::template patio_f<-N_two_pi>(2, 1);
			auto constexpr u = _op::diplo_f(+N_two)*_op::template patio_f<+N_two_pi>(2, 1);
			
			return S_::template method<Is...>(XTAL_REF_(o)*n)*u;
		};
		template <auto ...Is>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr n = _op::diplo_f(-N_two)*_op::template patio_f<-N_two_pi>(2, 1);
			auto constexpr u = _op::diplo_f(+N_two)*_op::template patio_f<+N_two_pi>(2, 1);
			
			return S_::template function<Is...>(XTAL_REF_(o)*n)*u;
		};

	};
};
template <>
struct dilated<0>
{
	template <class S>
	using subtype = S;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
