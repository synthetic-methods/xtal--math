#pragma once
#include "./any.hh"

#include "./root.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct near
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		/*/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		static_method(integral_variable_q auto n)
		noexcept -> XTAL_ALL_(n)
		{
			using U = XTAL_ALL_(n);
			using U_fix = bond::fixture<U>;

			--n;
			bond::seek_forward_f<bond::math::bit_ceiling_f(U_fix::full.depth)>([&]<constant_q I> (I) XTAL_0FN {n |= n >> (1 << I{});});
			++n;
			return n;
		}
		/***/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		static_method(real_variable_q auto u)
		noexcept -> XTAL_ALL_(u)
		{
			using U = XTAL_ALL_(u);
			using U_fix = bond::fixture<U>;
			using U_sigma = typename U_fix::sigma_type;
			using U_alpha = typename U_fix::alpha_type;

			U_alpha constexpr N_half_sqrt = root_f<2>(2.);
			U_sigma constexpr N_mask      = U_fix::sign.mask|U_fix::exponent.mask;
			return _xtd::bit_cast<U_alpha>(_xtd::bit_cast<U_sigma>(u*N_half_sqrt)&N_mask);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		static_method(complex_variable_q auto const &u)
		noexcept -> XTAL_ALL_(u)
		{
			return {static_method<Ns...>(u.real()), static_method<Ns...>(u.imag())};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
using    near_t = process::confined_t<near<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(return,inline,let)
near_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return near_t<Ms...>::static_method(XTAL_REF_(oo)...);
	return near_t<>::template static_method<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
