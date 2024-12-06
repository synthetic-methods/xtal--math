#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_two=0, int N_two_pi=0> struct   dilate;
template <int N_two=0, int N_two_pi=0> using    dilate_t = process::confined_t<dilate<N_two, N_two_pi>>;
template <int N_two=0, int N_two_pi=0>
XTAL_DEF_(short)
XTAL_LET dilate_f(auto &&o)
noexcept -> decltype(auto)
{
	return dilate_t<N_two, N_two_pi>::template function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int N_two, int N_two_pi>
struct dilate
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr n = _op::diplo_f(-N_two)*_op::template patio_f<-N_two_pi>(2, 1);
		//	auto constexpr u = _op::diplo_f(+N_two)*_op::template patio_f<+N_two_pi>(2, 1);
			
			return XTAL_REF_(o)*(n);
		};

	};
};
template <>
struct dilate<0, 0>
{
	template <class S>
	using subtype = S;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
