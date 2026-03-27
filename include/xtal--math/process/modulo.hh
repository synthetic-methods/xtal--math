#pragma once
#include "./any.hh"

#include "./part.hh"
#include "./nearest.hh"
#include "../atom/phason.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\returns The value wrapped within the range `{0, M}` if `M` is `unsigned`,
         or `{-2^M, +2^M}` if `M` is `signed`.
*/

template <auto ...Ms>
struct modulo;

template <auto ...Ms>
XTAL_TYP_(let) modulo_t = process::confined_t<modulo<Ms...>>;



////////////////////////////////////////////////////////////////////////////////

template <cardinal_q auto M_mod>
struct modulo<M_mod>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&z)
		noexcept -> auto
		requires un_v<atom::math::phason_q<decltype(z)>>
		{
			using bond::math::bit_count_f;
			using Z     = XTAL_ALL_(z);
			using Z_fit = bond::fit<Z>;
			XTAL_IF0
			XTAL_0IF (complex_field_q<Z>) {
				return complexion_f(method_f<Ns...>(z.real()), method_f<Ns...>(z.imag()));
			}
			XTAL_0IF (integral_constant_q<Z>) {
				return constant_t<method_f<Ns...>(z())>{};
			}
			XTAL_0IF (0 == bit_count_f(M_mod)) {
				return Z{};
			}
			XTAL_0IF (1 == bit_count_f(M_mod) and integral_q<Z>) {
				Z constexpr Z_mod{M_mod - 1};
				return XTAL_REF_(z)&Z_mod;
			}
			XTAL_0IF (2 <= bit_count_f(M_mod)) {
				Z constexpr Z_mod{M_mod};
				return (((XTAL_REF_(z)%Z_mod) + Z_mod)%Z_mod);
			}
		}

	};
};
template <ordinal_q auto M_mod>
struct modulo<M_mod>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&z)
		noexcept -> auto
		requires un_v<atom::math::phason_q<decltype(z)>>
		{
			using bond::math::bit_count_f;
			using Z     = XTAL_ALL_(z);
			using Z_fit = bond::fit<Z>;
			XTAL_IF0
			XTAL_0IF (complex_field_q<Z>) {
				return complexion_f(method_f<Ns...>(z.real()), method_f<Ns...>(z.imag()));
			}
			XTAL_0IF (integral_constant_q<Z>) {
				return constant_t<method_f<Ns...>(z())>{};
			}
			XTAL_0IF (1 <= M_mod and integral_q<Z>) {
				auto constexpr  Z_mod = Z_fit::sigma_1 << M_mod;
				return modulo_t<Z_mod>::template method_f<Ns...>(XTAL_REF_(z));
			}
			XTAL_0IF (M_mod <= 0 and integral_q<Z>) {
				return Z{};
			}
			XTAL_0IF (M_mod != 0) {
				auto constexpr Z_up = Z_fit::haplo_f(M_mod);
				auto constexpr Z_dn = Z_fit::diplo_f(M_mod);
				return modulo_t<-0>::template method_f<Ns...>(XTAL_REF_(z)*Z_up)*Z_dn;
			}
			XTAL_0IF_(to) (bond::math::bit_fraction_f(XTAL_REF_(z)))
			XTAL_0IF_(to) (z - nearest_f(z))
			XTAL_0IF_(void)
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&z)
		noexcept -> auto
		requires in_v<atom::math::phason_q<decltype(z)>>
		{
			using Z     = XTAL_ALL_(z);
			using Z_fit = bond::fit<Z>;
			XTAL_IF0
			XTAL_0IF (M_mod <  0) {
				auto constexpr N_mod = -M_mod;
				return (XTAL_REF_(z) << N_mod) >> N_mod;
			}
			XTAL_0IF (0 <= M_mod) {
				return XTAL_REF_(z);
			}
		}

	};
};
template <>
struct modulo<0, 0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_re=0, int N_im=0>
		XTAL_DEF_(return,set)
		method_f(auto &&x)
		noexcept -> decltype(auto)
		{
			using X = XTAL_ALL_(x);

			int constexpr K_re = N_re&1;
			int constexpr K_im = N_im&1;

			XTAL_IF0
			XTAL_0IF (atom::math::phason_complex_q<X> and K_re == 0 and K_im == 0) {
				auto const y = XTAL_REF_(x) (0);
				return _std::pair{y.real(), y.imag()};
			}
			XTAL_0IF (atom::math::phason_complex_q<X>) {
				using W = typename X::value_type;
				using V = typename W::value_type;
				V constexpr v = bond::fit<V>::sign.mask;
				return method_f(XTAL_REF_(x) - X{K_re*v, K_im*v});
			}
			XTAL_0IF (atom::math::complex_phason_q<X>) {
				using W = typename X::value_type;
				using V = typename W::value_type;
				V constexpr v = bond::fit<V>::sign.mask;
				auto x0_ = x.real(); if constexpr (K_re) x0_[0] -= v;
				auto xO_ = x.imag(); if constexpr (K_im) xO_[0] -= v;
				auto xO  = xO_(0);
				auto x0  = x0_(0);
				return _std::pair{XTAL_MOV_(x0), XTAL_MOV_(xO)};
			}
			XTAL_0IF (complex_variable_q<X>) {
				auto x0  = x.real(); if constexpr (K_re) x0 -= half*part_f<signed>(x0);
				auto xO  = x.imag(); if constexpr (K_im) xO -= half*part_f<signed>(xO);
				xO -= nearest_f<>(xO);
				x0 -= nearest_f<>(x0);
				return _std::pair{XTAL_MOV_(x0), XTAL_MOV_(xO)};
			}
			XTAL_0IF (simplex_variable_q<X> and K_re == 0) {
				return _std::pair{x - nearest_f(x), X{}};
			}
			XTAL_0IF (simplex_variable_q<X>) {
				auto const v = half*part_f<signed>(x);
				return method_f(XTAL_REF_(x) - XTAL_MOV_(v));
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
XTAL_DEF_(return,inline,let)
modulo_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return modulo_t<Ms...>::method_f(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
