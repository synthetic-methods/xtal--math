#pragma once
#include "./any.hh"

#include "./decompose.hh"
#include "./nearest.hh"
#include "../atom/phason.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\returns The argument modulo `0 < M_mod`,
or wrapped around `0` spanning `2^M_mod`.
*/

template <auto M_mod=null_type{}>
struct modulo;

template <auto M_mod=null_type{}>
XTAL_TYP_(let) modulo_t = process::confined_t<modulo<M_mod>>;



////////////////////////////////////////////////////////////////////////////////

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
				auto constexpr  Z_up  = Z_fit::haplo_f(M_mod);
				auto constexpr  Z_dn  = Z_fit::diplo_f(M_mod);
				auto constexpr  Z_mod = 0;
				return modulo_t<Z_mod>::template method_f<Ns...>(XTAL_REF_(z)*Z_up)*Z_dn;
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
template <fixed_shaped_q<null_type[2]> auto M_range>
struct modulo<M_range>
{
	static auto constexpr M_dn = get<0>(M_range);
	static auto constexpr M_up = get<1>(M_range);
	static auto constexpr M_delta = M_up - M_dn;
	static auto constexpr M_sigma = _xtd::make_unsigned_f(M_delta);
	static_assert(M_dn < M_up);

	using superkind = modulo<M_sigma>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&z)
		noexcept -> auto
		{
			using Z     = XTAL_ALL_(z);
			using Z_fit = bond::fit<Z>;

			XTAL_IF0
			XTAL_0IF (complex_field_q<Z>) {
				return complexion_f(method_f<Ns...>(z.real()), method_f<Ns...>(z.imag()));
			}
			XTAL_0IF_(else) {
				auto constexpr U_dn = Z_fit::alpha_f(M_dn);
				auto constexpr U_up = Z_fit::alpha_f(M_up);
				return S_::template method_f<Ns...>(XTAL_REF_(z) - U_dn) + U_dn;
			}
		};

	};
};


////////////////////////////////////////////////////////////////////////////////

template <>
XTAL_TYP_(new) modulo<null_type{}> : modulo<0> {};

template <auto M_mod=null_type{}, auto ...Ns>
XTAL_DEF_(let) modulo_f = [] XTAL_1FN_(call) (modulo_t<M_mod>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
