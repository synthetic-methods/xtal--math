#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "../../provision/prewarping.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Prepends a discrete/non-bandlimited control signal to the `method` arguments,

\tparam   M_end
Specifies the final `stage`.

When `M_end == 0`, an impulse is generated,
which is scaled by the sample-rate if `prewarping_q`.
This mode simulates a finite `DiracDelta[x]`.

When `M_end != 0`, a sustained pulse/gate is generated.
This mode simulates `HeavisidePi[# - 1/2]`.
*/
template <int M_end=1>
struct pulse
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
			auto const    &u_stage = S_::template head<stage_type>();
			auto const     i_stage = _xtd::make_unsigned_f(u_stage.head());
			auto constexpr I_stage = _xtd::make_unsigned_f(M_end);
			auto const     x_input = static_cast<typename state_type::value_type>(i_stage < I_stage);
			return S_::template method<Ns...>(x_input*half, XTAL_REF_(oo)...)*two;
		}

	};
};
template <>
struct pulse<0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;
		using T_ = typename S_::self_type;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::state_type;
		using typename S_::stage_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto u_warp, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto &u_stage = S_::template head<stage_type>();
			auto  x_stage = 0 == u_stage;
			auto  x_input = static_cast<typename state_type::value_type>(x_stage)*root_f<-1>(u_warp);
			static_assert(provision::math::prewarping_q<T_>);
			u_stage |= x_stage;
			return S_::template method<Ns...>(x_input*half, u_warp, XTAL_REF_(oo)...)*two;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
