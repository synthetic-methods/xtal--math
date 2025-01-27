#pragma once
#include "./any.hh"

#include "./unity.hh"
#include "../taylor/logarithm.hh"



XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines the pair `1^{#,-#} &`. \

///\note\
Pronounced "double-unity". \

template <int M_ism=0, typename ...As> requires in_n<M_ism, 0, 1, 2>
struct wnity
:	process::lift<wnity<M_ism>, bond::compose<As...>>
{
};
template <>
struct   wnity<>
{
	using limit_type = occur::math::limit_t<(1<<3)>;

	template <class S>
	using subtype = bond::compose_s<S, provision::context<void
	,	typename limit_type::template dispatch<>
	>>;

};
template <int M_ism=1, typename ...As>
using    wnity_t = process::confined_t<wnity<M_ism, As...>, wnity<>>;


////////////////////////////////////////////////////////////////////////////////

template <>
struct wnity<1> : wnity<>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(complex_field_q auto const &t)
		noexcept -> decltype(auto)
		{
			return static_method<N_lim>(t.real(), t.imag());
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(auto &&t_1, simplex_field_q auto &&t_i)
		noexcept -> decltype(auto)
		{
			auto constexpr exp = [] XTAL_0FN_(alias) (taylor::logarithm_t<-1, 1>::template static_method<2>);
			using T_i = XTAL_ALL_(t_i); using _fix = bond::fixture<T_i>;

			auto const o = static_method<N_lim>(XTAL_REF_(t_1));
			auto const e = exp(XTAL_REF_(t_i)*_fix::patio_f(-2));
			return o*roots_t<1>::static_method(e);
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(simplex_field_q auto &&t_1)
		noexcept -> decltype(auto)
		{
			auto const o = objective_f(unity_t<1>::template static_method<N_lim>(XTAL_REF_(t_1)));
			auto const p = complexion_f(o.real(),  o.imag());
			auto const q = complexion_f(o.real(), -o.imag());
			//\
			return atom::couple_f(o, conj(o));
			return atom::couple_f(p, q);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
