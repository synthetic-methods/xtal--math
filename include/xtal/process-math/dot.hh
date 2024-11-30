#pragma once
#include "./any.hh"

#include "./root.hh"




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
noexcept -> decltype(auto)
{
	return root_f<M_pow, M_zap>(norm(XTAL_REF_(o)));
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
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&o)
		noexcept -> auto
		{
			return dot_f<M_pow, M_zap>(XTAL_REF_(o));
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
