#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Provides result-expansion (w.r.t. argument-reduction) for `unity`, \
i.e. squaring for both circular and hyperbolic results.

template <int M_ism=1> struct   disunity;
template <int M_ism=1> using    disunity_t = process::confined_t<disunity<M_ism>>;

template <int M_ism=1, auto ...Ns>
XTAL_DEF_(return,inline,let)
disunity_f(auto &&o)
noexcept -> decltype(auto)
{
	return disunity_t<M_ism>::template static_method<Ns...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 0>
struct disunity<M_ism>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		static_method(auto &&o)
		noexcept -> decltype(auto)
		{
			return XTAL_REF_(o);
		}

	};
};
template <int M_ism> requires in_n<M_ism, 1, 2>
struct disunity<M_ism>
{
	static constexpr int I_sgn = sign_n<(M_ism&1)^1, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		static_method(complex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using U    = XTAL_ALL_(o);
			using U_fix = bond::fixture<U>;

			auto constexpr i_sgn = U_fix::alpha_f(I_sgn);

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
