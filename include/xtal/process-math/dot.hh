#pragma once
#include "./any.hh"

#include "./root.hh"
#include "./square.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_pow=1, int M_zap=-1> requires in_n<M_pow, 1, 2,-1,-2>
XTAL_TYP dot;

template <int M_pow=1, int M_zap=-1> requires in_n<M_pow, 1, 2,-1,-2>
XTAL_USE dot_t = process::confined_t<dot<M_pow, M_zap>>;

template <int M_pow=1, int M_zap=-1> requires in_n<M_pow, 1, 2,-1,-2>
XTAL_DEF_(return,inline)
XTAL_LET dot_f(auto &&o)
XTAL_0EX -> decltype(auto)
{
	return dot_t<M_pow, M_zap>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_pow, int M_zap> requires in_n<M_pow, 1, 2,-1,-2>
struct dot//<M_pow, M_zap>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_SET function(simplex_field_q auto &&o)
		XTAL_0EX -> auto
		{
			return root_f<M_pow, M_zap>(square_f(XTAL_REF_(o)));
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_SET function(complex_field_q auto &&o)
		XTAL_0EX -> auto
		{
			return root_f<M_pow, M_zap>(square_f(o.real(), o.imag()));
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
