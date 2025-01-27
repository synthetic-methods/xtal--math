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
XTAL_DEF_(return,inline,let)
roots_f(auto &&o)
noexcept -> decltype(auto)
{
	return roots_t<M_exp, M_cut>::static_method(XTAL_REF_(o));
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
		XTAL_DEF_(return,inline,set)
		static_method(auto &&w)
		noexcept -> auto
		{
			using    _fix = bond::fixture<XTAL_ALL_(w)>;
			auto const o =   objective_f(XTAL_REF_(w));

			XTAL_IF0
			XTAL_0IF (M_exp ==  1) {auto const q = root_f<-1, M_cut>(o); return atom::couple_f(o,       q);}
			XTAL_0IF (M_exp ==  2) {auto const q = root_f<-2, M_cut>(o); return atom::couple_f(o*q,     q);}
			XTAL_0IF (M_exp ==  3) {auto const q = root_f<-3, M_cut>(o); return atom::couple_f(o*q*q,   q);}
			XTAL_0IF (M_exp ==  4) {auto const q = root_f<-4, M_cut>(o); return atom::couple_f(o*q*q*q, q);}
			XTAL_0IF (M_exp == -1) {auto const q = root_f<-1, M_cut>(o); return atom::couple_f(q,       o);}
			XTAL_0IF (M_exp == -2) {auto const q = root_f<-2, M_cut>(o); return atom::couple_f(q,     q*o);}
			XTAL_0IF (M_exp == -3) {auto const q = root_f<-3, M_cut>(o); return atom::couple_f(q,   q*q*o);}
			XTAL_0IF (M_exp == -4) {auto const q = root_f<-4, M_cut>(o); return atom::couple_f(q, q*q*q*o);}
			XTAL_0IF (0 <  M_exp) {auto const q = root_f<-M_exp, M_cut, Ns...>(o); return atom::couple_f(o *power_f<+M_exp - 1>(q), q);}
			XTAL_0IF (M_exp <  0) {auto const q = root_f<+M_exp, M_cut, Ns...>(o); return atom::couple_f(q, power_f<-M_exp - 1>(q)* o);}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
