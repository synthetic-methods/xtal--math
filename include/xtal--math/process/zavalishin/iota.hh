#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "./prewarped.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Prepends a discrete/non-bandlimited control signal to the `method` arguments, \
with the final stage specified by `M_end`. \

///\
When `M_end == 0`, an impulse/trigger is generated, \
which is scaled by the sample-rate if `prewarped_q`. \
This mode simulates a finite `DiracDelta[x]`. \

///\
When `M_end != 0`, a gate/hold is generated. \
This mode simulates `HeavisidePi[# - 1/2]`. \

template <int M_end=1>
struct iota
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::stage_type;
		using typename S_::input_type;

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
			auto const     x_input = static_cast<input_type>(i_stage < I_stage);
			return S_::template method<Ns...>(x_input, XTAL_REF_(oo)...);
		}

	};
};
template <>
struct iota<0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;
		using T_ = typename S_::self_type;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::stage_type;
		using typename S_::input_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto s_scale, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto &u_stage = S_::template head<stage_type>();
			auto  x_stage = 0 == u_stage;
			auto  x_input = static_cast<input_type>(x_stage);
			if constexpr (prewarped_q<T_>) {
				x_input *= root_f<-1>(s_scale);
			}
			u_stage |= x_stage;
			return S_::template method<Ns...>(x_input, s_scale, XTAL_REF_(oo)...);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
