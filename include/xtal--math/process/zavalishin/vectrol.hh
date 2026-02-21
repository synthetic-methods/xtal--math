#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "../taylor/logarithm.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Scales `damp` by the dot-product of the internal state and supplied `reshape`.
\note    Input is restricted to `U_pole` because the filter-state is managed out-of-band.
\todo    Perform `dot` with state-derived outputs, defining the high-pass as `one - s_shape_.sum()`.
\todo    Include input-gain?
*/
template <auto ...As>	struct  vectrol;
template <auto ...As>	using   vectrol_t = confined_t<vectrol<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <auto ...As>
struct vectrol
{
	using archetype = occur::context_t<vectrol>;
	using superkind = typename archetype::template attach<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;

		XTAL_DEF_(set) exp_f = [] XTAL_1FN_(call) (taylor::logarithm_t<-1>::template method_f<3>);

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::shape_parameter;
		using typename S_::shape_type;
		using typename S_::state_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto const &[s_state_] = S_::template memory<state_type>();
			auto const & s_shape_  = S_::template head<shape_parameter>().head();
			return S_::template method<Ns...>(XTAL_REF_(oo)...,
				exp_f(-dot_f(s_state_, s_shape_)));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <auto ..._s>
struct context<process::math::zavalishin::vectrol<_s...>>
{
	using superkind = context<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
	
	public:
		using S_::S_;

		template <extent_type N_mask=1>
		struct   attach
		{
			template <class R>
			using subtype = bond::compose_s<R, typename S_::template   attach<N_mask>
			,	provision::voiced<void
				,	typename R::shape_parameter::template   attach<N_mask>
				>
			>;

		};

	};
};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
