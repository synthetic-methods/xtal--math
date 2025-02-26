#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "../taylor/logarithm.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Scales `damping` by the dot-product of the internal state and supplied `reshape`. \

///\note\
Input is restricted to `U_pole` because the filter-state is managed out-of-band. \

template <typename ...As>	struct  vectrol;
template <typename ...As>	using   vectrol_t = process::confined_t<vectrol<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct any<vectrol<As...>> : any<filter<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

template <vector_q A>
struct vectrol<A>
{
	using metakind = any<vectrol<A>>;

	using   state_type = typename metakind::   state_type;
	using   curve_type = typename metakind::   curve_type;
	using recurve_type = typename metakind:: recurve_type;
	
	using superkind = bond::compose<bond::tag<vectrol_t>
	,	typename recurve_type::template attach<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto x_input, auto s_scale, auto s_damping, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto constexpr abs = [] XTAL_1FN_(call) (taylor::logarithm_t<-1>::template method_f<0>);

			auto const [s_]       = S_::template memory<  state_type>();
			auto const &s_reshape = S_::template   head<recurve_type>().head();
			auto const  s_product = dot_f(s_, s_reshape);
			
			s_damping *= abs(s_product*half);

			return S_::template method<Ns...>(x_input, s_scale, s_damping, XTAL_REF_(oo)...);
		}

	};
};
template <scalar_q A>
struct vectrol<A>
:	vectrol<A[2]>
{
};
template <>
struct vectrol<>
:	vectrol<typename bond::fit<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
