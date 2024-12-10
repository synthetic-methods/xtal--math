#pragma once
#include "./any.hh"
#include "./root.hh"





XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_exp=1, int M_cut=0> struct   roots;
template <int M_exp=1, int M_cut=0> using    roots_t = process::confined_t<roots<M_exp>>;
template <int M_exp=1, int M_cut=0>
XTAL_DEF_(short)
XTAL_LET roots_f(auto &&o)
noexcept -> decltype(auto)
{
	return roots_t<M_exp, M_cut>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_exp, int M_cut>
struct roots
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&w)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(w)>;
			auto constexpr _1 = _op::alpha_1;

			auto const o = objective_f(XTAL_REF_(w));

			XTAL_IF0
			XTAL_0IF (M_exp ==  1) {auto const q = root_f<-1, M_cut>(o); return duple_f(o,       q);}
			XTAL_0IF (M_exp ==  2) {auto const q = root_f<-2, M_cut>(o); return duple_f(o*q,     q);}
			XTAL_0IF (M_exp ==  3) {auto const q = root_f<-3, M_cut>(o); return duple_f(o*q*q,   q);}
			XTAL_0IF (M_exp ==  4) {auto const q = root_f<-4, M_cut>(o); return duple_f(o*q*q*q, q);}
			XTAL_0IF (M_exp == -1) {auto const q = root_f<-1, M_cut>(o); return duple_f(q,       o);}
			XTAL_0IF (M_exp == -2) {auto const q = root_f<-2, M_cut>(o); return duple_f(q,     q*o);}
			XTAL_0IF (M_exp == -3) {auto const q = root_f<-3, M_cut>(o); return duple_f(q,   q*q*o);}
			XTAL_0IF (M_exp == -4) {auto const q = root_f<-4, M_cut>(o); return duple_f(q, q*q*q*o);}
			XTAL_0IF (0 <  M_exp) {auto const q = root_f<-M_exp, M_cut, Ns...>(o); return duple_f(o *_op::template explo_f<+M_exp - 1>(q), q);}
			XTAL_0IF (M_exp <  0) {auto const q = root_f<+M_exp, M_cut, Ns...>(o); return duple_f(q, _op::template explo_f<-M_exp - 1>(q)* o);}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
