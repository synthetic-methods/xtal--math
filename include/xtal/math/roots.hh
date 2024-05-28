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
	XTAL_LET_(int) K_pow = -bond::operate<int>::designed_f(N_pow);

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
			using _std::sqrt;

			using op = bond::operate<decltype(o)>;
			auto constexpr _1 = op::alpha_1;
			/*/
			auto const q = op::template root_f<K_pow, N_lim>(o);
			XTAL_IF0
			XTAL_0IF (N_pow ==  2) {return bond::couple_f(o*q, q);}
			XTAL_0IF (N_pow ==  1) {return bond::couple_f(o,   q);}
			XTAL_0IF (N_pow == -1) {return bond::couple_f(q,   o);}
			XTAL_0IF (N_pow == -2) {return bond::couple_f(q, q*o);}
			/*/
			XTAL_IF0
			XTAL_0IF (N_pow ==  2) {auto const q = _1/sqrt(o); return bond::couple_f(o*q, q);}
			XTAL_0IF (N_pow ==  1) {auto const q = _1/    (o); return bond::couple_f(o,   q);}
			XTAL_0IF (N_pow == -1) {auto const q = _1/    (o); return bond::couple_f(q,   o);}
			XTAL_0IF (N_pow == -2) {auto const q = _1/sqrt(o); return bond::couple_f(q, q*o);}
			/***/
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
