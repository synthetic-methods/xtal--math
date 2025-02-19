#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "./gate.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Manages a ringing filter, \
scaling `damping` according to the stored `stage`. \

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
	using superkind = bond::tag<ring_t>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&x_input, auto s_scale, auto s_damping, auto &&...oo)
		noexcept -> decltype(auto)
		{
			using S_stage = typename S_::stage_type;
			using U = XTAL_ALL_(s_damping);
			U constexpr N0{0};
			U constexpr N1{1};
			U constexpr N2{2};
			U constexpr N2_up = root_f< 2>(N2);
			U constexpr N2_dn = root_f<-2>(N2);

			_std::array<U, 4> constexpr stage_damping{N0, N2_dn, N1, N2_up};
			s_damping *= stage_damping[0b11&S_::template head<S_stage>()];

			return S_::template method<Ns...>(XTAL_REF_(x_input), s_scale, s_damping, XTAL_REF_(oo)...);
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
