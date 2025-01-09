#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1> struct   squishy;
template <int M_ism=1> using    squishy_t = process::confined_t<squishy<M_ism>>;

template <int M_ism=1, auto ...Ns>
XTAL_DEF_(short)
XTAL_LET squishy_f(auto &&o)
noexcept -> decltype(auto)
{
	return squishy_t<M_ism>::template function<Ns...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 0>
struct squishy<M_ism>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			return XTAL_REF_(o);
		}

	};
};
template <int M_ism> requires in_n<M_ism, 1, 2>
struct squishy<M_ism>
{
	static constexpr int I_sgn = signum_n<(M_ism&1)^1, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using U    = XTAL_ALL_(o);
			using U_op = bond::operate<U>;

			XTAL_LET i_sgn = U_op::alpha_f(I_sgn);

			auto y = o.imag();
			auto x = o.real();
			auto const xx = x*x;
			auto const yy = y*y;
			auto const y_ = two*x*y;
			auto const x_ = term_f(xx, yy, +i_sgn);
			auto const w_ = term_f(xx, yy, -i_sgn);
			auto const m_ = root_f<-1>(w_);
			return complexion_f(x_, y_)*(m_);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
