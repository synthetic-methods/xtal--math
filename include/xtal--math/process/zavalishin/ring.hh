#pragma once
#include "./any.hh"

#include "./filter.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Manages a ringing filter with stored `damping` and `balance`.

///\note\
Input is restricted to `U_pole` because the filter-state is managed out-of-band. \

template <typename ...As> struct   ring;
template <typename ...As> using    ring_t = process::confined_t<ring<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class U_pole, int N_pole>
struct ring<U_pole[N_pole]>
:	filter<U_pole[N_pole]>
{
	using control_type = absolve_u<U_pole>;
	using damping_type = occur::inferred_t<struct DAMPING, control_type>;
	using balance_type = occur::inferred_t<struct BALANCE, control_type>;

	using superkind = bond::compose<bond::tag<ring_t>
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
		XTAL_FX4_(alias) (XTAL_DEF_(return,inline,let) damping(auto &&...oo), S_::template head<damping_type>(XTAL_REF_(oo)...))
		XTAL_FX4_(alias) (XTAL_DEF_(return,inline,let) balance(auto &&...oo), S_::template head<balance_type>(XTAL_REF_(oo)...))

	public:// FLUX

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion> requires in_n<N_ion, +1>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			auto constexpr N_dn = root_f<2>(control_type{two});
			auto constexpr N_up = power_f<-2*7>(N_dn);

			if constexpr (occur::stage_q<decltype(o)>) {
				switch (o.head()) {
					case  0: {(void) damping(-N_up); break;}
					case  1: {(void) damping( N_up); break;}
					case -1: {(void) damping( N_dn); break;}
				}
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion> requires in_n<N_ion, -1>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			using _fix = bond::fixture<control_type>;
			if constexpr (occur::stage_q<decltype(o)>) {
				switch (o.head()) {
					case  0: {break;}
					case  1: {break;}
					case -1: {
						auto constexpr N_thresh = _fix::haplo_f(7);
						auto constexpr M_thresh = square_f(N_thresh);
						using          U_poles_ = atom::couple_t<U_pole[N_pole]>;
						//\note\
						Because `N_ord` may be less than `N_pole`, \
						access to `U_poles_` requres resetting `stow` when `N_ord` changes. \

						//\note\
						May need to finess checking the zeroness of `states_`, \
						since the components are scaled by powers of gain. \

						auto const [states_] = S_::template stow<U_poles_>();
						return bond::pack_dot_f(states_) < M_thresh and S_::template fuse<N_ion>(XTAL_REF_(o));
					}
				}
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}

	public:
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(U_pole s_scale, auto &&...oo)
		noexcept -> decltype(auto)
		{
			U_pole x_input{one};
			if constexpr XTAL_TRY_(do) (x_input *= S_::sample().rate())
			auto const s_damping = static_cast<U_pole>(S_::template head<damping_type>());
			auto const y_balance = static_cast<U_pole>(S_::template head<balance_type>());
			return S_::template method<Ns...>(U_pole{one}, s_scale, s_damping, y_balance) - term_f<-1, 2>(one, y_balance)*(x_input);
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
:	ring<typename bond::fixture<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
