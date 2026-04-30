#pragma once
#include "./any.hh"

#include "../../process/dot.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Manages the lifecycle of the current voice.
\todo    Provision toggle to select between monophonic/polyphonic.
*/
template <auto  ..._s>	XTAL_TYP_(new) reuse   : bond::compose<reuse<_s>...> {};
template <auto  ..._s>	XTAL_TYP_(set) reuse_t = confined_t<reuse<_s...>>;
template <class ..._s>	XTAL_TYP_(ask) reuse_q = bond::tab_inner_p<reuse<>, _s...>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Responds to `efflux(occur::stage_f(+0))` by clearing the filter state.
*/
template <>
struct reuse<>
{
};
template <int M_ind>
struct reuse<M_ind>
{
	using superkind = bond::tab<reuse<>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FLOW

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		/*!
		\brief   Resets the filter state on `stage`-onset.
		\note    Active only when `M_ind == -0`.
		*/
		template <signed N_ion> requires in_v<M_ind, -0> and in_v<N_ion, -1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			if (o.head() == 0) {
				S_::memory(constant_t<>{});
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		/*!
		\brief   Reports the `state` of the filter via `stage` event, with `-1` indicating "off".
		\note    Active only when `M_ind == -1`.
		\todo    Threshold is currently fixed at `2^-16`, or around `-48dBFS`, but should be configurable.
		\returns `1` if the state is under threshold, `0` otherwise.
		*/
		template <signed N_ion> requires in_v<M_ind, -1> and in_v<N_ion, +1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			using U_order =  typename S_     ::order_attribute;
			using U_state =  typename S_     :: data_type;
			using V_scale =  typename U_state::scale_type;
			using V_fit   = bond::fit<V_scale>;
			auto constexpr esquilon = V_fit::haplo_f(0x11);// ~-48dBFS squared...

			signed x = S_::template fuse<N_ion>(XTAL_REF_(o));
			auto const [states_] = S_::template memory<U_state>();
			auto const n_order   = U_order{S_::self()};
			if (o == -1 and 1 <= n_order) {
				V_scale sum{};
				for (int i{}; i < n_order; ++i) {
					sum += process::math::dot_f(states_[i]);
				}
				x &= sum < esquilon;
			}
			return x;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
