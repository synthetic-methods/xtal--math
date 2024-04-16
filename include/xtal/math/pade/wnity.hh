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

template <int M_ism=1, typename ...As> struct wnity {static_assert(M_ism);};
template <int M_ism=1, typename ...As> using  wnity_t = process::confined_t<wnity<M_ism, As...>>;

template <int M_ism, bond::compose_q ...As> requires some_q<As...>
struct wnity<M_ism, As...>: process::chain<wnity<M_ism>, As...> {};


////////////////////////////////////////////////////////////////////////////////

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
		XTAL_FN2 function(auto &&t_re, simplex_field_q auto &&t_im)
		XTAL_0EX
		{
		//	using T_re = XTAL_TYP_(t_re);
			using T_im = XTAL_TYP_(t_im);
			auto constexpr f = bond::realize<T_im>::patio_f(-2);
			auto const     o = unity_t<1>::template function<N_lim>(XTAL_REF_(t_re));
			return algebra::scalar_f(o, _std::conj(o))*
				roots_t<1>::template function(_std::exp(f*XTAL_REF_(t_im)));
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
