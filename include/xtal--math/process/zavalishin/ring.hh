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
struct ring<U_pole[N_pole]> : filter<U_pole[N_pole]>
{
	using U_fit = bond::fit<U_pole>;

	using filter_kind = filter<U_pole[N_pole]>;
	using  order_type = typename filter_kind::order_type;

	using control_type = absolve_u<U_pole>;
	using damping_type = occur::inferred_t<struct DAMPING, control_type>;
	using balance_type = occur::inferred_t<struct BALANCE, control_type>;

	using stage_type = occur::stage_t<>;

	using superkind = bond::compose<bond::tag<ring_t>
	,	typename damping_type::template attach<>
	,	typename balance_type::template attach<>
	,	typename   stage_type::template expect<>
	,	filter<U_pole[N_pole]>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using U_poles_ = atom::couple_t<U_pole[N_pole]>;

	public:// CONSTRUCT
		using S_::S_;

	public:
	//	TODO: Need parareter mapping/restriction on creation/update...
		XTAL_FX4_(to) (XTAL_DEF_(return,inline,let) damping(auto &&...oo), S_::template head<damping_type>(XTAL_REF_(oo)...))
		XTAL_FX4_(to) (XTAL_DEF_(return,inline,let) balance(auto &&...oo), S_::template head<balance_type>(XTAL_REF_(oo)...))
		XTAL_FX4_(to) (XTAL_DEF_(return,inline,let)   stage(auto &&...oo), S_::template head<  stage_type>(XTAL_REF_(oo)...))

	public:// FLUX

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
			auto const [states_] = S_::template memory<U_poles_>();
			signed x = S_::template fuse<N_ion>(XTAL_REF_(o));

			switch (o.head()) {
			//\
			case  0: (void) damping(          (U_fit::minilon_f(0))); break;
			case  0: (void) damping(          (U_fit::  alpha_f(0))); break;
			case  1: (void) damping(root_f<-2>(U_fit::  haplo_f(1))); break;
			case -1: (void) damping(root_f< 2>(U_fit::  diplo_f(1))); break;
			}
			return x;
		}
		template <signed N_ion> requires in_n<N_ion, -1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			auto const [states_] = S_::template memory<U_poles_>();
			signed x = S_::template fuse<N_ion>(XTAL_REF_(o));
			
			switch (o.head()) {
			case  0: break;
			case  1: break;
			case -1:
			//	TODO: Accommodate frequency scaling when calculating the product?

				auto constexpr N_thresh_sqr = bond::fit<control_type>::haplo_f(7 + 1);
				unsigned int const order = S_::template head<order_type>().head();
				U_pole y{};
				for (unsigned int i{}; i < order; ++i) {
					y = term_f<1, 2>(y, states_[i]);
				}
				x &= y < N_thresh_sqr;
			}
			return x;
		}

	public:
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(U_pole s_scale, auto &&...oo)
		noexcept -> decltype(auto)
		{
			//\
			auto const x_input   = control_type{one};
			auto const x_input   = static_cast<U_pole>(0 <= stage());
			auto const s_damping = static_cast<U_pole>(S_::template head<damping_type>());
			auto const y_balance = static_cast<U_pole>(S_::template head<balance_type>());
			return S_::template method<Ns...>(x_input, s_scale, s_damping, y_balance) - term_f<-1, 2>(one, y_balance)*(x_input);
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
