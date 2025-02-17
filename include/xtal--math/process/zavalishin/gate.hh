#pragma once
#include "./any.hh"

#include "./filter.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Prepends an input signal to `method`'s arguments: \
`1` while the current `stage >= 0`, \
`0` otherwise. \

template <typename ...As>	struct  gate;
template <typename ...As>	using   gate_t = process::confined_t<gate<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <vector_q A>
struct gate<A>
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

	public:// ACCESS
		using S_::self;

	public:// OPERATE
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto const &u_stage = S_::template head<stage_type>();
			auto const  x_input = static_cast<valued_u<state_type>>(0 <= u_stage);
			return S_::template method<Ns...>(x_input, XTAL_REF_(oo)...);
		}

	};
};
template <scalar_q A>
struct gate<A>
:	gate<A[2]>
{
};
template <>
struct gate<>
:	gate<typename bond::fit<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
