#pragma once
#include "./any.hh"

#include "./root.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_exp=1, int M_cut=0> requires in_q<M_exp, 1, 2,-1,-2>
struct   dot;

template <int M_exp=1, int M_cut=0> requires in_q<M_exp, 1, 2,-1,-2>
using    dot_t = process::confined_t<dot<M_exp, M_cut>>;

template <int M_exp=1, int M_cut=0> requires in_q<M_exp, 1, 2,-1,-2>
XTAL_DEF_(short)
XTAL_LET dot_f(auto &&o)
noexcept -> decltype(auto)
{
	return root_f<M_exp, M_cut>(norm(XTAL_REF_(o)));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_exp, int M_cut> requires in_q<M_exp, 1, 2,-1,-2>
struct dot//<M_exp, M_cut>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> auto
		{
			return dot_f<M_exp, M_cut>(XTAL_REF_(o));
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
