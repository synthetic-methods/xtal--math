#pragma once
#include "./any.hh"

#include "./filter.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Manages a ringing filter, \
scaling `damping` according to the stored `stage`. \

///\note\
Input is restricted to `U_pole` because the filter-state is managed out-of-band. \

////////////////////////////////////////////////////////////////////////////////

template <int M_end=0>
struct ring
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;

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
			s_damping *= _std::array<U, 4>{N0, N2_dn, N1, N2_up}[0b11&S_::template head<S_stage>()];

			return S_::template method<Ns...>(XTAL_REF_(x_input), s_scale, s_damping, XTAL_REF_(oo)...);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
