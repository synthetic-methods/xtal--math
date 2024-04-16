#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_pow=1> struct roots;
template <int N_pow=1> using  roots_t = process::confined_t<roots<N_pow>>;
template <int N_pow=1>
XTAL_FN2 roots_f(auto &&o)
XTAL_0EX
{
	return roots_t<N_pow>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int N_pow>
struct roots
{
	XTAL_LET_(int) K_pow = -bond::realized::designed_f(N_pow);

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>// requires sign_p<N_pow, 0>
		XTAL_FN2 function(auto const &o)
		XTAL_0EX
		{
			using re = bond::realize<decltype(o)>;
			auto const n = re::template root_f<K_pow, N_lim>(o);
			XTAL_IF0
			XTAL_0IF_(N_pow ==  2) {return algebra::scalar_f(o*n, n);}
			XTAL_0IF_(N_pow ==  1) {return algebra::scalar_f(o,   n);}
			XTAL_0IF_(N_pow == -1) {return algebra::scalar_f(n,   o);}
			XTAL_0IF_(N_pow == -2) {return algebra::scalar_f(n, n*o);}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
