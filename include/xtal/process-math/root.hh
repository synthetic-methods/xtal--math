#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_pow=1, int M_zap=-1> requires inclusive_q<M_pow, 1, 2,-1,-2>
XTAL_TYP root;

template <int M_pow=1, int M_zap=-1> requires inclusive_q<M_pow, 1, 2,-1,-2>
XTAL_USE root_t = process::confined_t<root<M_pow, M_zap>>;

template <int M_pow=1, int M_zap=-1> requires inclusive_q<M_pow, 1, 2,-1,-2>
XTAL_DEF_(return,inline)
XTAL_LET root_f(auto &&o)
XTAL_0EX -> decltype(auto)
{
	return root_t<M_pow, M_zap>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_pow, int M_zap> requires inclusive_q<M_pow, 1, 2,-1,-2>
struct root//<M_pow, M_zap>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&o)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<XTAL_ALL_(o)>;
			XTAL_IF0
			XTAL_0IF (M_pow ==  1) {return                   _op::template punctured_f<M_zap*1>(XTAL_REF_(o)) ;}
			XTAL_0IF (M_pow ==  2) {return              sqrt(_op::template punctured_f<M_zap*2>(XTAL_REF_(o)));}
			XTAL_0IF (M_pow == -1) {return _op::alpha_1/    (_op::template punctured_f<M_zap*1>(XTAL_REF_(o)));}
			XTAL_0IF (M_pow == -2) {return _op::alpha_1/sqrt(_op::template punctured_f<M_zap*2>(XTAL_REF_(o)));}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
