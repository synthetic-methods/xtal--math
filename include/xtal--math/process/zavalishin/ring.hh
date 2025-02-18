#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "./gate.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Manages a ringing filter with stored `damping`, \
which is dynamically configured by the `influx`ed `stage`. \

///\note\
Input is restricted to `U_pole` because the filter-state is managed out-of-band. \

template <class ...As>	struct  ring;
template <class ...As>	using   ring_t = process::confined_t<ring<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct any<ring<As...>> : any<filter<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

template <vector_q A>
struct ring<A>
{
	using _fit = bond::fit<A>;

	using     metakind = any<ring<A>>;
	using   state_type = typename metakind::   state_type;
	using damping_type = typename metakind:: damping_type;

	using superkind = bond::compose<bond::tag<ring_t>
	,	typename damping_type::template attend<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
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
		template <signed N_ion> requires in_n<N_ion,  1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			auto const [states_] = S_::template memory<state_type>();
			signed x = S_::template fuse<N_ion>(XTAL_REF_(o));

			switch (o.head()) {
			case  0: (void) S_::template head<damping_type>(          (_fit::alpha_0)); break;
			case  1: (void) S_::template head<damping_type>(root_f< 2>(_fit::haplo_1)); break;
			case -1: (void) S_::template head<damping_type>(root_f<-2>(_fit::haplo_1)); break;
			}
			return x;
		}

	};
};
template <scalar_q A>
struct ring<A>
:	ring<A[2]>
{
};
template <>
struct ring<>
:	ring<typename bond::fit<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
