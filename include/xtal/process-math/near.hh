#pragma once
#include "./any.hh"






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

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_number_q auto u)
		noexcept -> XTAL_ALL_(u)
		{
			using U = XTAL_ALL_(u);
			using U_op = bond::operate<U>;
			using U_sigma = typename U_op::sigma_type;
			using U_alpha = typename U_op::alpha_type;

			XTAL_LET N_half = U_op::alpha_f(1.4142135623730950488016887242096980786L);// Sqrt[2]
			XTAL_LET N_mask = U_op::sign.mask|U_op::exponent.mask;
			return _xtd::bit_cast<U_alpha>(_xtd::bit_cast<U_sigma>(u*N_half)&N_mask);
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_number_q auto const &u)
		noexcept -> XTAL_ALL_(u)
		{
			return {function<Ns...>(u.real()), function<Ns...>(u.imag())};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
using    near_t = process::confined_t<near<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET near_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return near_t<Ms...>::function(XTAL_REF_(oo)...);
	return near_t<>::template function<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
