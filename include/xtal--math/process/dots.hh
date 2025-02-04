#pragma once
#include "./any.hh"

#include "./roots.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_exp=1, int M_cut=0> requires in_n<M_exp, 1, 2,-1,-2>
struct   dots;

template <int M_exp=1, int M_cut=0> requires in_n<M_exp, 1, 2,-1,-2>
using    dots_t = process::confined_t<dots<M_exp, M_cut>>;

template <int M_exp=1, int M_cut=0> requires in_n<M_exp, 1, 2,-1,-2>
XTAL_DEF_(return,inline,let)
dots_f(auto &&o)
noexcept -> decltype(auto)
{
	return roots_f<M_exp, M_cut>(norm(XTAL_REF_(o)));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_exp, int M_cut> requires in_n<M_exp, 1, 2,-1,-2>
struct dots//<M_exp, M_cut>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			return dots_f<M_exp, M_cut>(XTAL_REF_(o));
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
