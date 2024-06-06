#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_two=0, int N_two_pi=0>
struct dilating;


////////////////////////////////////////////////////////////////////////////////

template <int N_two, int N_two_pi>
struct dilating
{
	using subkind = bond::tag<process::link>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto &&o)
		XTAL_0EX
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr n = _op::diplo_f(-N_two)*_op::template patio_f<-N_two_pi>(2, 1);
			auto constexpr u = _op::diplo_f(+N_two)*_op::template patio_f<+N_two_pi>(2, 1);
			
			return S_::template function<Is...>(XTAL_REF_(o)*(n))*(u);
		};

	};
};
template <>
struct dilating<0>
{
	template <class S>
	using subtype = S;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
