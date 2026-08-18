#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Evaluates `Total[{##}^2] &` (using fused multiply-add, if supported by the compiler).
*/
template <auto ...Ms>	struct  square;
template <auto ...Ms>	using   square_t = process::confined_t<square<Ms...>>;

template <int N_alt=1, int N_sig=1, int N_sqr=1>
XTAL_VAL_(return,inline,let)
square_f(auto const &x, auto &&...xs)
noexcept -> auto
{
	using X = XTAL_ALL_(x);

	XTAL_IF0
	XTAL_0IF (2 <= N_sqr) {
		return square_f<N_alt, N_sig, N_sqr - 1>(square_f<N_alt, N_sig>(x, XTAL_REF_(xs)...));
	}
	XTAL_0IF (1 == N_sqr) {
		XTAL_IF0
		XTAL_0IF (1 <= sizeof...(xs)) {
			using W = unstruct_t<X>;
			auto constexpr K_sig = N_alt*N_sig;
			auto constexpr k_sig =   (W) K_sig;
			return xtd::plus_multiplies{} (k_sig*square_f<N_alt, K_sig>(XTAL_REF_(xs)...), x, x);
		}
		XTAL_0IF (complex_field_q<X>) {
			using W = unstruct_t<X>;
			auto constexpr K_sig = N_alt*N_sig;
			auto constexpr k_sig =   (W) K_sig;
		//	auto const &[x_re, x_im] = destruct_f(x);
			auto const x_re = x.real();
			auto const x_im = x.imag();
			return complexion_f(square_f<N_alt, -1>(x_re, x_im), (W)two*x_re*x_im);
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
		XTAL_VAL_(return,inline,let)
		method(auto &&...xs) const
		noexcept -> decltype(auto)
		{
			return square_f<Ns...>(XTAL_REF_(xs)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
