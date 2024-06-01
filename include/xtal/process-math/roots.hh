#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
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
	XTAL_LET_(int) K_pow = -bond::operate<int>::designed_f(N_pow);

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>// requires sign_p<N_pow, 0>
		XTAL_FN2 function(auto &&w)
		XTAL_0EX
		{
			using _std::sqrt;

			using Op = bond::operate<decltype(w)>;
			auto constexpr _1 = Op::alpha_1;

			auto const o = objective_f(XTAL_REF_(w));

			XTAL_IF0
			XTAL_0IF (N_pow ==  2) {auto const q = _1/sqrt(o); return algebra::scalar_f(o*q, q);}
			XTAL_0IF (N_pow ==  1) {auto const q = _1/    (o); return algebra::scalar_f(o,   q);}
			XTAL_0IF (N_pow == -1) {auto const q = _1/    (o); return algebra::scalar_f(q,   o);}
			XTAL_0IF (N_pow == -2) {auto const q = _1/sqrt(o); return algebra::scalar_f(q, q*o);}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
