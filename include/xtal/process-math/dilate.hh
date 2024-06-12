#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_two=0, int N_two_pi=0> XTAL_TYP dilate;
template <int N_two=0, int N_two_pi=0> XTAL_USE dilate_t = process::confined_t<dilate<N_two, N_two_pi>>;
template <int N_two=0, int N_two_pi=0>
XTAL_LET dilate_f = [] (auto &&o)
XTAL_0FN {
	using _op = bond::operate<decltype(o)>;
	auto constexpr n = _op::diplo_f(-N_two)*_op::template patio_f<-N_two_pi>(2, 1);
//	auto constexpr u = _op::diplo_f(+N_two)*_op::template patio_f<+N_two_pi>(2, 1);
	
	return XTAL_REF_(o)*(n);
};

////////////////////////////////////////////////////////////////////////////////

template <int N_two, int N_two_pi>
struct dilate
:	process::defer<decltype(dilate_f<N_two, N_two_pi>)>
{
	
//	template <class S>
//	class subtype: public bond::compose_s<S>
//	{
//		using S_ = bond::compose_s<S>;
//
//	public:
//		using S_::S_;
//
//		template <auto ...>
//		XTAL_DEF_(return,inline,static)
//		XTAL_RET function(auto &&o)
//		XTAL_0EX
//		{
//			using _op = bond::operate<decltype(o)>;
//			auto constexpr n = _op::diplo_f(-N_two)*_op::template patio_f<-N_two_pi>(2, 1);
//		//	auto constexpr u = _op::diplo_f(+N_two)*_op::template patio_f<+N_two_pi>(2, 1);
//			
//			return XTAL_REF_(o)*(n);
//		};
//
//	};
};
template <>
struct dilate<0>
{
	template <class S>
	using subtype = S;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
