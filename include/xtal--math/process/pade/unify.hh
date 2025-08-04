#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides result-expansion (w.r.t. argument-reduction) for `unity`.

Implements squaring for both circular and hyperbolic results.
*/
template <int M_ism=1>
struct unify;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 0>
struct unify<M_ism>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return XTAL_REF_(o);
		}

	};
};
template <int M_ism> requires in_n<M_ism, 1, 2>
struct unify<M_ism>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using U    = XTAL_ALL_(o);
			using U_fit = bond::fit<U>;

			auto constexpr N_sgn = U_fit::alpha_f(cosign_v<M_ism>);

			auto y = o.imag();
			auto x = o.real();
			auto const xx = x*x;
			auto const yy = y*y;
			auto const y_ = two*x*y;
			auto const x_ = term_f(xx, yy, +N_sgn);
			auto const w_ = term_f(xx, yy, -N_sgn);
			auto const m_ = root_f<-1>(w_);
			return complexion_f(x_, y_)*(m_);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1>
using unify_t = process::confined_t<unify<M_ism>>;

template <int M_ism=1, int ...Ns>
XTAL_DEF_(let)
unify_f = [] XTAL_1FN_(call) (unify_t<M_ism>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
