#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_two=0, int N_two_pi=0> struct dilate;
template <int N_two=0, int N_two_pi=0> using  dilate_t = process::confined_t<dilate<N_two, N_two_pi>>;


////////////////////////////////////////////////////////////////////////////////

template <int N_two, int N_two_pi>
struct dilate
{
	using subkind = bond::tag<process::chain>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using re = bond::realize<decltype(o)>;
			auto constexpr n = re::diplo_f(-N_two)*re::template patio_f<-N_two_pi>(2, 1);
			auto constexpr u = re::diplo_f(+N_two)*re::template patio_f<+N_two_pi>(2, 1);
			
			return S_::template function<Is...>(XTAL_REF_(o)*(n))*(u);
		};

	};
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
