#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "./trigger.hh"
#include "../taylor/logarithm.hh"


XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Manages a negative-ramping filter with stored `damping`.

///\note\
Input is restricted to `U_pole` because the filter-state is managed out-of-band. \

template <typename ...As>	struct  roll;
template <typename ...As>	using   roll_t = process::confined_t<roll<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct any<roll<As...>> : any<filter<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

template <vector_q A>
struct roll<A>
{
	using _fit = bond::fit<A>;

	using     metakind = any<roll<A>>;
	using   order_type = typename metakind::   order_type;
	using   state_type = typename metakind::   state_type;
	using   curve_type = typename metakind::   curve_type;
	using recurve_type = typename metakind:: recurve_type;
	using damping_type = typename metakind:: damping_type;
	
	using superkind = bond::compose<bond::tag<roll_t>
	,	typename recurve_type::template attach<>
	,	typename damping_type::template attach<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::state_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto x_input, auto s_scale, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto constexpr warp_f = [] XTAL_1FN_(function) (taylor::logarithm_t<-1>::template method_f<0>);

			auto const [s_]       = S_::template memory<  state_type>();
			auto const &s_recurve = S_::template   head<recurve_type>().head();
			auto const  s_damping = S_::template   head<damping_type>()*warp_f(dot_f(s_, s_recurve));

			if constexpr (prewarped_q<T_> and trigger_q<T_>) {
				x_input *= root_f<-1>(s_scale);
			}
			return S_::template method<Ns...>(x_input, s_scale, s_damping, XTAL_REF_(oo)...);
		}

	};
};
template <scalar_q A>
struct roll<A>
:	roll<A[2]>
{
};
template <>
struct roll<>
:	roll<typename bond::fit<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
