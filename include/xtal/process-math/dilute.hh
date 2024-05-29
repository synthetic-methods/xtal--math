#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_two=0, int N_two_pi=0> struct dilute;
template <int N_two=0, int N_two_pi=0> using  dilute_t = process::confined_t<dilute<N_two, N_two_pi>>;


////////////////////////////////////////////////////////////////////////////////

template <int N_two, int N_two_pi>
struct dilute
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using op = bond::operate<decltype(o)>;
			auto constexpr n = op::diplo_f(-N_two)*op::template patio_f<-N_two_pi>(2, 1);
		//	auto constexpr u = op::diplo_f(+N_two)*op::template patio_f<+N_two_pi>(2, 1);
			
			return XTAL_REF_(o)*(n);
		};

	};
};
template <>
struct dilute<0>
{
	template <class S>
	using subtype = S;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
