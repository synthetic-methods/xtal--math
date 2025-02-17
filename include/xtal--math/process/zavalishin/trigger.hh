#pragma once
#include "./any.hh"

#include "./filter.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Prepends an input signal to `method`'s arguments: \
`1` for the first sample after receiving `occur::stage_f(0)`, \
`0` otherwise. \

///\note\
Advances the current stage to `1` once the first sample has been processed. \

template <typename ...As>	struct  trigger;
template <typename ...As>	using   trigger_t = process::confined_t<trigger<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <vector_q A>
struct trigger<A>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::stage_type;
		using typename S_::state_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto      &u_stage = S_::template head<stage_type>();
			auto const x_stage = 0 == u_stage;
			auto const x_input = static_cast<valued_u<state_type>>(x_stage);
			u_stage |= x_stage;
			return S_::template method<Ns...>(x_input, XTAL_REF_(oo)...);
		}

	};
};
template <scalar_q A>
struct trigger<A>
:	trigger<A[2]>
{
};
template <>
struct trigger<>
:	trigger<typename bond::fit<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
