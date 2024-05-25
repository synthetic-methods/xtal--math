#pragma once
#include "./any.hh"
#include "./unity.hh"
#include "../roots.hh"




XTAL_ENV_(push)
namespace xtal::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the pair `1^{#,-#} &`. \

///\note\
Pronounced "double-unity". \

template <int M_ism=1, typename ...As> XTAL_TYP wnity {static_assert(M_ism);};
template <int M_ism=1, typename ...As> XTAL_USE wnity_t = process::confined_t<wnity<M_ism, As...>>;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, bond::compose_q ...As> requires some_q<As...>
struct wnity<M_ism, As...>
:	process::chain<wnity<M_ism>, As...>
{
};
template <>
struct wnity<1>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_FN2 function(complex_field_q auto const &t)
		XTAL_0EX
		{
			return function<N_lim>(t.real(), t.imag());
		}
		template <int N_lim=-1>
		XTAL_FN2 function(auto &&t_1, simplex_field_q auto &&t_i)
		XTAL_0EX
		{
			using _std::exp;

			using T_i = XTAL_TYP_(t_i); using op = bond::operate<T_i>;

			return function(XTAL_REF_(t_1))*
				roots_t<1>::template function(exp(XTAL_REF_(t_i)*op::patio_f(-2)));
		}
		template <int N_lim=-1>
		XTAL_FN2 function(simplex_field_q auto &&t_1)
		XTAL_0EX
		{
			using _std::conj;

			auto const o = unity_t<1>::template function<N_lim>(XTAL_REF_(t_1));
			return bond::couple_f(o, conj(o));
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
