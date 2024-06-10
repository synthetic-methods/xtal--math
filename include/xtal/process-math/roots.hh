#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_pow=1> XTAL_TYP roots;
template <int N_pow=1> XTAL_USE roots_t = process::confined_t<roots<N_pow>>;
template <int N_pow=1>
XTAL_DEF_(return,inline)
XTAL_REF roots_f(auto &&o)
XTAL_0EX
{
	return roots_t<N_pow>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int N_pow>
struct roots
{
	static constexpr int K_pow = -bond::operate<int>::designed_f(N_pow);

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_REF function(auto &&w)
		XTAL_0EX
		{
			using _std::sqrt;

			using _op = bond::operate<decltype(w)>;
			auto constexpr _1 = _op::alpha_1;

			auto const o = objective_f(XTAL_REF_(w));

			XTAL_IF0
			XTAL_0IF (N_pow ==  2) {auto const q = _1/sqrt(o); return duple_f(o*q, q);}
			XTAL_0IF (N_pow ==  1) {auto const q = _1/    (o); return duple_f(o,   q);}
			XTAL_0IF (N_pow == -1) {auto const q = _1/    (o); return duple_f(q,   o);}
			XTAL_0IF (N_pow == -2) {auto const q = _1/sqrt(o); return duple_f(q, q*o);}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
