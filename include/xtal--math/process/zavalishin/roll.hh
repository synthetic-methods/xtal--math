#pragma once
#include "./any.hh"

#include "./filter.hh"




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
	using   state_type = typename metakind::   state_type;
	using   style_type = typename metakind::   style_type;
	using restyle_type = typename metakind:: restyle_type;
	using damping_type = typename metakind:: damping_type;
	
	using superkind = bond::compose<bond::tag<roll_t>
	,	typename restyle_type::template attach<>
	,	typename damping_type::template attach<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::state_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto x_input, auto s_scale, auto &&...oo)
		noexcept -> decltype(auto)
		{
		//	NOTE: Assumes `s_scale` is prewarped...
		//	s_scale *= S_::sampling().period();
			x_input *= root_f<-1>(s_scale);// Assuming `prewarped`...

		//	TODO: Abstract `warp_f` as `shaper_f`...
			auto constexpr warp_f = [] (auto &&x) XTAL_0FN_(to) (root_f<2>(term_f<1, 2>(one, x)) + x);

			auto const & f_      = S_::template   head<restyle_type>();
			auto const & t_      = f_ .template   head<  style_type>();
			auto const  [s_]     = S_::template memory<  state_type>();
			auto const y_damping = S_::template   head<damping_type>()*warp_f(dot_f(s_, t_));

			return S_::template method<Ns...>(x_input, s_scale, y_damping, XTAL_REF_(oo)...);
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
