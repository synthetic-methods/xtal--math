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
struct  wnity<>
{
	using limit_type = occur::inferred_t<union LIMIT, bond::seek_s<(1<<3)>>;

	template <class S>
	using subtype = bond::compose_s<S, provision::voiced<void
	,	typename limit_type::template dispatch<>
	>>;

};
template <int M_ism=1, typename ...As>
using   wnity_t = process::confined_t<wnity<M_ism, As...>, wnity<>>;


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
			static_assert(real_variable_q<decltype(t_re)>);
			static_assert(real_variable_q<decltype(t_im)>);

			auto constexpr exp = [] XTAL_1FN_(call) (taylor::logarithm_t<-1, 1>::template method_f<2>);
			using _fit = bond::fit<decltype(t_re), decltype(t_im)>;

			return method_f<N_lim>(XTAL_REF_(t_re))*
				roots_t<1>::method_f(exp(XTAL_REF_(t_im)*_fit::patio_f(-2)));
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&t_re)
		noexcept -> decltype(auto)
		{
			auto const p = objective_f(unity_t<1>::template method_f<N_lim>(XTAL_REF_(t_re)));
			auto const q = objective_f(conj(p));
			return atom::couple_f(p, q);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
