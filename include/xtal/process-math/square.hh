#pragma once
#include "./any.hh"
#include "./term.hh"





XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Evaluates `Total[{x, xs}^2]` (using fused multiply-add, if supported by the compiler). \

template <auto ...Ms> XTAL_TYP square;
template <auto ...Ms> XTAL_USE square_t = process::confined_t<square<Ms...>>;

template <int N_alt=1, int N_sgn=1>
XTAL_DEF_(return,inline)
XTAL_LET square_f(auto const &x, auto &&...xs)
noexcept -> auto
{
	using X = XTAL_ALL_(x);
	using _op = bond::operate<X>;
	
	auto constexpr K_sgn = N_alt*N_sgn;

	XTAL_IF0
	XTAL_0IF (1 <= sizeof...(xs)) {
		return term_f(_op::alpha_f(K_sgn)*square_f<N_alt, K_sgn>(XTAL_REF_(xs)...), x, x);
	}
	XTAL_0IF (complex_field_q<X>) {
	//	auto const &[x_re, x_im] = involved_f(x);
		auto const x_re = x.real();
		auto const x_im = x.imag();
		return complexion_f(square_f<-N_alt>(x_re, x_im), _op::diplo_1*x_re*x_im);
	}
//	XTAL_0IF (simplex_field_q<X>) {
	XTAL_0IF_(else) {
		return x*x;
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
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&...xs)
		noexcept -> decltype(auto)
		{
			return square_f<Ns...>(XTAL_REF_(xs)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
