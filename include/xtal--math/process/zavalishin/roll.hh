#pragma once
#include "./any.hh"

#include "./filter.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Manages a negative-ramping filter with stored `damping` and `balance`.

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
	using   scope_type = typename metakind::   scope_type;
	using rescope_type = typename metakind:: rescope_type;
	using damping_type = typename metakind:: damping_type;
	using balance_type = typename metakind:: balance_type;
	
	using superkind = bond::compose<bond::tag<roll_t>
	,	typename rescope_type::template attach<>
	,	typename damping_type::template attach<>
	,	typename balance_type::template attend<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::state_type;

	public:// ACCESS
	//	TODO: Need parareter mapping/restriction on creation/update...

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,let)
		damping(auto &&...oo), S_::template head<damping_type>(XTAL_REF_(oo)...))

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,let)
		balance(auto &&...oo), S_::template head<balance_type>(XTAL_REF_(oo)...))

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

			auto const & f_      = S_::template   head<rescope_type>();
			auto const & t_      = f_ .template   head<  scope_type>();
			auto const  [s_]     = S_::template memory<  state_type>();
			auto const &[s0, s1] = s_;
			auto const y_damping = warp_f(s1 + dot_f<-1>(s_, t_))*damping();

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
