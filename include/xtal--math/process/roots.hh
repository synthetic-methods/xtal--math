#pragma once
#include "./any.hh"
#include "./root.hh"





XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_exp=1, int M_cut=0>	struct  roots;
template <int M_exp=1, int M_cut=0>	using   roots_t = process::confined_t<roots<M_exp>>;
template <int M_exp=1, int M_cut=0>
XTAL_DEF_(return,inline,let)
roots_f(auto &&o)
noexcept -> decltype(auto)
{
	return roots_t<M_exp, M_cut>::method_f(XTAL_REF_(o));
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
		method_f(auto &&w)
		noexcept -> auto
		requires different_q<decltype(w), decltype(objective_f(w))>
		{
			return method_f(objective_f(XTAL_REF_(w)));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&w)
		noexcept -> auto
		requires      same_q<decltype(w), decltype(objective_f(w))>
		{
			using W0 = XTAL_ALL_(w);
			using W2 = atom::couple_t<W0[2]>;
			XTAL_IF0
			XTAL_0IF (M_exp ==  1) {auto const n = root_f<-1, M_cut>(w); return W2{w,       n};}
			XTAL_0IF (M_exp ==  2) {auto const n = root_f<-2, M_cut>(w); return W2{w*n,     n};}
			XTAL_0IF (M_exp ==  3) {auto const n = root_f<-3, M_cut>(w); return W2{w*n*n,   n};}
			XTAL_0IF (M_exp ==  4) {auto const n = root_f<-4, M_cut>(w); return W2{w*n*n*n, n};}
			XTAL_0IF (M_exp == -1) {auto const n = root_f<-1, M_cut>(w); return W2{n,       w};}
			XTAL_0IF (M_exp == -2) {auto const n = root_f<-2, M_cut>(w); return W2{n,     n*w};}
			XTAL_0IF (M_exp == -3) {auto const n = root_f<-3, M_cut>(w); return W2{n,   n*n*w};}
			XTAL_0IF (M_exp == -4) {auto const n = root_f<-4, M_cut>(w); return W2{n, n*n*n*w};}
			XTAL_0IF (0 <   M_exp) {auto const n = root_f<-M_exp, M_cut, Ns...>(w); return W2{w*power_f<+M_exp - 1>(n),n};}
			XTAL_0IF (M_exp <   0) {auto const n = root_f<+M_exp, M_cut, Ns...>(w); return W2{n,power_f<-M_exp - 1>(n)*w};}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
