#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Manages the lifecycle of the current voice.
\todo    Incorporate toggle to select between monophonic/polyphonic.
\todo    Consider making reset/reuse a passive `provision`,
         and move this into `meta` (or something).
*/
template <auto  ...Ms>	struct  reuse;
template <auto  ...Ms>	using   reuse_t = process::confined_t<reuse<Ms...>>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Responds to `efflux(occur::stage_f(+0))` by clearing the filter state.
*/
template <>
struct reuse< 0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS
		using S_::self;

	public:// FLOW

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion> requires in_v<N_ion, -1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			if (o.head() == 0) {
				S_::memory(constant_t<>{});
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Responds to `influx(occur::stage_f(-1))` by returning `1` if the state is under threshold, `0` otherwise.
\node    Threshold is currently fixed at `2^-16`, or around `-48dBFS`.
\todo    Allow configurable threshold.
\todo    Accommodate gain when calculating the `dot` with state?
*/
template <>
struct reuse<-1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::order_attribute;
		using typename S_::state_type;

	public:// ACCESS
		using S_::self;

	public:// FLOW

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(flow::assessing_q<occur::stage_t<>> auto &&o)
		noexcept -> signed
		{
			return fuse<+1>(XTAL_REF_(o).tail());
		}
		template <signed N_ion> requires in_v<N_ion, +1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			auto const [states_] = S_::template memory<state_type>();
			using  Y = XTAL_ALL_(dot_f(states_));
			signed x = S_::template fuse<N_ion>(XTAL_REF_(o));

			auto const ord = order_attribute{self()};
			if (o.head() == -1 and 1 <= ord) {
				Y dot{};
				for (int i{}; i < ord; ++i) {
					dot = term_f(XTAL_MOV_(dot), dot_f(states_[i]));
				}
				x &= dot < bond::fit<Y>::haplo_f(0x11);// Epsilon (~-48dBFS) squared...
			}
			return x;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
