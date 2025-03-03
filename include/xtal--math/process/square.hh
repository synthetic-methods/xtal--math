#pragma once
#include "./any.hh"
#include "./term.hh"





XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Evaluates `Total[{##}^2] &` (using fused multiply-add, if supported by the compiler). \

template <auto ...Ms>	struct  square;
template <auto ...Ms>	using   square_t = process::confined_t<square<Ms...>>;

template <int N_sqr=1, int N_alt=1, int N_sgn=1>
XTAL_DEF_(return,inline,let)
square_f(auto const &x, auto &&...xs)
noexcept -> auto
{
	XTAL_IF0
	XTAL_0IF (2 <= N_sqr) {
		return square_f<N_sqr - 1, N_alt, N_sgn>(square_f<1, N_alt, N_sgn>(x, XTAL_REF_(xs)...));
	}
	XTAL_0IF (1 == N_sqr) {
		using X    = XTAL_ALL_(x);
		using X_fit = bond::fit<X>;
		
		auto constexpr K_sgn = N_alt*N_sgn;

		XTAL_IF0
		XTAL_0IF (1 <= sizeof...(xs)) {
			return term_f(X_fit::alpha_f(K_sgn)*square_f<1, N_alt, K_sgn>(XTAL_REF_(xs)...), x, x);
		}
		XTAL_0IF (complex_field_q<X>) {
		//	auto const &[x_re, x_im] = destruct_f(x);
			auto const x_re = x.real();
			auto const x_im = x.imag();
			return complexion_f(square_f<1, -N_alt>(x_re, x_im), X_fit::diplo_1*x_re*x_im);
		}
	//	XTAL_0IF (simplex_field_q<X>) {
		XTAL_0IF_(else) {
			return x*x;
		}
	}
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct square
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&...xs)
		noexcept -> decltype(auto)
		{
			return square_f<Ns...>(XTAL_REF_(xs)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
