#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
struct near
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_number_q auto u)
		noexcept -> XTAL_ALL_(u)
		{
			using U = XTAL_ALL_(u);
			using U_op = bond::operate<U>;
			auto const N_shift = U_op::exponential_f(u);// TODO: Optimize...
			auto const U_up = U_op::haplo_f(N_shift);
			auto const U_dn = U_op::diplo_f(N_shift);

			u *= U_up;

			XTAL_IF0
			XTAL_0IF_(consteval) {
				u += U_op::haplo_1;
				u = static_cast<U>(static_cast<integer_type>(u));
			}
			XTAL_0IF_(else) {
				u = round(u);
			}
			u *= U_dn;
			return u;
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_number_q auto const &u)
		noexcept -> XTAL_ALL_(u)
		{
			return {function<Ns...>(u.real()), function<Ns...>(u.imag())};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
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
