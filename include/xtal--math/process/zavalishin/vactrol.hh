#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "../pade/tangy.hh"
#include "../taylor/logarithm.hh"


XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Scales the frequency using the input/state difference and the supplied `reshape`. \

///\note\
Input is restricted to `U_pole` because the filter-state is managed out-of-band. \

template <typename ...As>	struct  vactrol;
template <typename ...As>	using   vactrol_t = process::confined_t<vactrol<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct any<vactrol<As...>> : any<filter<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

template <vector_q A>
struct vactrol<A>
{
	using metakind = any<vactrol<A>>;

	using   state_type = typename metakind::   state_type;
	using   shape_type = typename metakind::   shape_type;
	using reshape_type = typename metakind:: reshape_type;
	
	using superkind = bond::compose<bond::tag<vactrol_t>
	,	typename reshape_type::template attach<>
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
		method(auto x_input, auto s_scale, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto constexpr abs = [] XTAL_1FN_(call) (taylor::logarithm_t<-1>::template method_f<0>);

			auto const &[d0, d1] = S_::template head<reshape_type>().head();
			auto const  [s_]     = S_::template memory<state_type>();

			auto const v = pade::tangy_t<1>::template method_f< 1>(half*d1);
			auto const w = abs(x_input - s_.sum());
			s_scale /= d0*term_f(v, one - v, w);

			return S_::template method<Ns...>(x_input, s_scale, XTAL_REF_(oo)...);
		}

	};
};
template <scalar_q A>
struct vactrol<A>
:	vactrol<A[2]>
{
};
template <>
struct vactrol<>
:	vactrol<typename bond::fit<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
