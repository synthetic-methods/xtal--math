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
		XTAL_DEF_(short,static)
		XTAL_LET function(integral_variable_q auto n)
		noexcept -> XTAL_ALL_(n)
		{
			using U = XTAL_ALL_(n);
			using U_op = bond::operate<U>;

			--n;
			bond::seek_forward_f<U_op::bit_floor_f(U_op::full.depth)>([&]<constant_q I> (I) XTAL_0FN {n |= n >> (1 << I{});});
			++n;
			return n;
		}
		/***/
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_variable_q auto u)
		noexcept -> XTAL_ALL_(u)
		{
			using U = XTAL_ALL_(u);
			using U_op = bond::operate<U>;
			using U_sigma = typename U_op::sigma_type;
			using U_alpha = typename U_op::alpha_type;

			U_alpha constexpr N_half = root_f<2>(2.);
			U_sigma constexpr N_mask = U_op::sign.mask|U_op::exponent.mask;
			return _xtd::bit_cast<U_alpha>(_xtd::bit_cast<U_sigma>(u*N_half)&N_mask);
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_variable_q auto const &u)
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
