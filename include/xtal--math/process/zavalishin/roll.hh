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
Manages a negative-ramping filter, \
scaling `damping` by the product of `recurve` and the internal state.

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
	using metakind = any<roll<A>>;

	using   state_type = typename metakind::   state_type;
	using   curve_type = typename metakind::   curve_type;
	using recurve_type = typename metakind:: recurve_type;
	
	using superkind = bond::compose<bond::tag<roll_t>
	,	typename recurve_type::template attach<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto x_input, auto s_scale, auto s_damping, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto const [s_]       = S_::template memory<  state_type>();
			auto const &s_recurve = S_::template   head<recurve_type>();
			auto const  s_product = dot_f(s_, s_recurve.head());
			
			s_damping *= taylor::logarithm_t<-1>::template method_f<0>(s_product);

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
