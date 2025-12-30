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

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::shape_parameter;
		using typename S_::shape_type;
		using typename S_::state_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto x, auto s_gain, auto s_damp, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto const  [s_]     = S_::template memory<state_type     >();
			auto const & s_shape = S_::template   head<shape_parameter>();
			s_damp *= taylor::logarithm_t<-1>::template method_f<0>(dot_f(s_, s_shape.head()));
			return S_::template method<Ns...>(x, s_gain, s_damp, XTAL_REF_(oo)...);
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
				//\
				,	typename T_::shape_parameter::template   attach<N_mask>
				,	typename R ::shape_parameter::template   attach<N_mask>
				>
			>;

		};

	};
};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
