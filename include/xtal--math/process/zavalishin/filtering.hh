#pragma once
#include "./any.hh"

#include "../../provision/shaper.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Manages a ringing filter with stored `damping` and `balance`.

///\note\
Input is restricted to `U_pole` because the filter-state is managed out-of-band. \

template <typename ...As> struct   filtering;
template <typename ...As> using    filtering_t = process::confined_t<filtering<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class U_pole, int N_pole>
struct filtering<U_pole[N_pole]>
{
	using control_type = absolve_u<U_pole>;
	using damping_type = occur::inferred_t<struct DAMPING, control_type>;
	using balance_type = occur::inferred_t<struct BALANCE, control_type>;

	using superkind = bond::compose<bond::tag<filtering_t>
	,	typename damping_type::template attach<>
	,	typename balance_type::template attach<>
	,	filter<U_pole[N_pole]>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:
		XTAL_TO4_(XTAL_DEF_(let) damping(auto &&...oo), S_::template head<damping_type>(XTAL_REF_(oo)...))
		XTAL_TO4_(XTAL_DEF_(let) balance(auto &&...oo), S_::template head<balance_type>(XTAL_REF_(oo)...))

	public:// FLUX
		XTAL_DEF_(short)
		XTAL_LET infuse(auto &&o)
		noexcept -> signed
		{
			XTAL_LET N_dn = root_f<2>(control_type{two});
			XTAL_LET N_up = power_f<-2*7>(N_dn);

			if constexpr (occur::stage_q<decltype(o)>) {
				switch (o.head()) {
					case  0: {(void) damping(-N_up); break;}
					case  1: {(void) damping( N_up); break;}
					case -1: {(void) damping( N_dn); break;}
				}
			}
			return S_::infuse(XTAL_REF_(o));
		}
		XTAL_DEF_(short)
		XTAL_LET effuse(auto &&o)
		noexcept -> signed
		{
			using _op = bond::operate<control_type>;
			if constexpr (occur::stage_q<decltype(o)>) {
				switch (o.head()) {
					case  0: {break;}
					case  1: {break;}
					case -1: {
						XTAL_LET N_thresh = _op::haplo_f(7);
						XTAL_LET M_thresh = square_f(N_thresh);
						using    U_poles_ = arrange::couple_t<U_pole[N_pole]>;
						//\note\
						Because `N_ord` may be less than `N_pole`, \
						access to `U_poles_` requres resetting the `cache` when `N_ord` changes. \

						//\note\
						May need to finess checking the zeroness of `states_`, \
						since the components are scaled by powers of gain. \

						auto const [states_] = S_::template cache<U_poles_>();
						return bond::pack_dot_f(states_) < M_thresh and S_::effuse(XTAL_REF_(o));
					}
				}
			}
			return S_::effuse(XTAL_REF_(o));
		}

	public:
		template <auto ...Ns>
		XTAL_DEF_(short)
		XTAL_LET method(auto &&x_input, real_q auto s_scale)
		noexcept -> decltype(auto)
		{
			static_assert(same_q<U_pole, decltype(x_input)>);
		//	static_assert(N_top == 0);// Necessary?

			using X_op = bond::operate<decltype(x_input)>;
			auto const s_damping = X_op::alpha_f(S_::template head<damping_type>());
			auto const y_balance = X_op::alpha_f(S_::template head<balance_type>());
			return S::template method<Ns...>(XTAL_REF_(x_input), s_scale, s_damping, y_balance);
		}

	};
};
template <>
struct filtering<>
:	filtering<typename bond::operating::aphex_type[4]>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
