#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides the nearest power-of-two.
*/
template <auto ...Ms>
struct near;


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
		XTAL_DEF_(return,inline,let)
		method(integral_variable_q auto n)
		const noexcept -> XTAL_ALL_(n)
		{
			using U = XTAL_ALL_(n);
			using U_fit = bond::fit<U>;

			--n;
			bond::seek_to_e<bond::math::bit_ceiling_f(U_fit::full.depth)>([&]<constant_q I> (I) XTAL_0FN {n |= n >> (1 << I{});});
			++n;
			return n;
		}
		/***/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(real_variable_q auto u)
		const noexcept -> XTAL_ALL_(u)
		{
			using U = XTAL_ALL_(u);
			using U_fit = bond::fit<U>;
			using U_sigma = typename U_fit::sigma_type;
			using U_alpha = typename U_fit::alpha_type;

			U_alpha constexpr N_half_sqrt = _std::numbers::sqrt2_v<U_alpha>;// root_f<2>(2.);
			U_sigma constexpr N_mask      = U_fit::sign.mask|U_fit::exponent.mask;
			return _xtd::bit_cast<U_alpha>(_xtd::bit_cast<U_sigma>(u*N_half_sqrt)&N_mask);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(complex_variable_q auto const &u)
		const noexcept -> XTAL_ALL_(u)
		{
			return {method<Ns...>(u.real()), method<Ns...>(u.imag())};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
using   near_t = process::confined_t<near<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(return,inline,let)
near_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return near_t<Ms...>{}.method(XTAL_REF_(oo)...);
	return near_t<>{}.template method<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
