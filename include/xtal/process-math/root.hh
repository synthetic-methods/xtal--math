#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_pow=1> XTAL_TYP root;
template <int M_pow=1> XTAL_USE root_t = process::confined_t<root<M_pow>>;
template <int M_pow=1>
XTAL_DEF_(return,inline)
XTAL_REF root_f(auto &&o)
XTAL_0EX
{
	return root_t<M_pow>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_pow>
struct root//<M_pow>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_REF function(auto &&o)
		XTAL_0EX
		{
			using _op = bond::operate<XTAL_TYP_(o)>;
			XTAL_IF0
			XTAL_0IF (M_pow ==  2) {return              sqrt(XTAL_REF_(o));}
			XTAL_0IF (M_pow ==  1) {return     XTAL_TYP_(o) (XTAL_REF_(o));}
			XTAL_0IF (M_pow == -1) {return _op::alpha_1/    (XTAL_REF_(o));}
			XTAL_0IF (M_pow == -2) {return _op::alpha_1/sqrt(XTAL_REF_(o));}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
