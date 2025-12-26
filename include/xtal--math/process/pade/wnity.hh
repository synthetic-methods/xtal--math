#pragma once
#include "./any.hh"

#include "./unity.hh"
#include "../taylor/logarithm.hh"
#include "../taylor/octarithm.hh"


XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines the pair `1^{#,-#} &`.
\note    Pronounced "double-unity".
*/
template <int M_ism=1, int M_car=0>
struct wnity;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, int M_car> requires in_v<M_ism, 1> and in_v<M_car, 0, 1>
struct wnity<M_ism, M_car>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto const &t)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(t.real(), t.imag());
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&t_re, simplex_field_q auto &&t_im)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(XTAL_REF_(t_re))*
				roots_t<1>::method_f(taylor::octarithm_f<-1>(XTAL_REF_(t_im)));
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&t_re)
		noexcept -> decltype(auto)
		{
			auto const p = objective_f(unity_t<M_ism, M_car>::template method_f<N_lim>(XTAL_REF_(t_re)));
			auto const q = objective_f(conj(p));
			return atom::couple_f(p, q);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=1>
using wnity_t = process::confined_t<wnity<M_ism, M_car>>;

template <int M_ism=1, int M_car=1, int N_lim=2>
XTAL_DEF_(let)
wnity_f = [] XTAL_1FN_(call) (wnity_t<M_ism, M_car>::template method_f<N_lim>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
