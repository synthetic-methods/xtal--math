#pragma once
#include "./any.hh"

#include "./exponential.hh"
#include "./monomial.hh"
#include "./nearest.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct prime;

template <auto ...Ms>
using prime_t = process::confined_t<prime<Ms...>>;

template <auto ...Ms>
XTAL_DEF_(return,inline,let)
prime_f(auto &&...xs)
noexcept -> decltype(auto)
{
	return prime_t<Ms...>::method(XTAL_REF_(xs)...);
};


namespace _detail
{///////////////////////////////////////////////////////////////////////////////

//	Table[BitShiftRight[Prime[i + 1]] - i, {i, 0, 127}]
unsigned char constexpr U8_PRIME_HALF_SUBDEX[]
{	0x01U, 0x00U, 0x00U, 0x00U, 0x01U, 0x01U, 0x02U, 0x02U
,	0x03U, 0x05U, 0x05U, 0x07U, 0x08U, 0x08U, 0x09U, 0x0BU
,	0x0DU, 0x0DU, 0x0FU, 0x10U, 0x10U, 0x12U, 0x13U, 0x15U
,	0x18U, 0x19U, 0x19U, 0x1AU, 0x1AU, 0x1BU, 0x21U, 0x22U
,	0x24U, 0x24U, 0x28U, 0x28U, 0x2AU, 0x2CU, 0x2DU, 0x2FU
,	0x31U, 0x31U, 0x35U, 0x35U, 0x36U, 0x36U, 0x3BU, 0x40U
,	0x41U, 0x41U, 0x42U, 0x44U, 0x44U, 0x48U, 0x4AU, 0x4CU
,	0x4EU, 0x4EU, 0x50U, 0x51U, 0x51U, 0x55U, 0x5BU, 0x5CU
,	0x5CU, 0x5DU, 0x63U, 0x65U, 0x69U, 0x69U, 0x6AU, 0x6CU
,	0x6FU, 0x71U, 0x73U, 0x74U, 0x76U, 0x79U, 0x7AU, 0x7DU
,	0x81U, 0x81U, 0x85U, 0x85U, 0x87U, 0x88U, 0x8AU, 0x8DU
,	0x8EU, 0x8EU, 0x8FU, 0x94U, 0x97U, 0x98U, 0x9BU, 0x9CU
,	0x9EU, 0xA3U, 0xA3U, 0xABU, 0xADU, 0xB1U, 0xB3U, 0xB5U
,	0xB5U, 0xB7U, 0xBBU, 0xBDU, 0xBFU, 0xBFU, 0xC1U, 0xC3U
,	0xC4U, 0xC4U, 0xC9U, 0xCDU, 0xCDU, 0xCEU, 0xD0U, 0xD2U
,	0xD2U, 0xD7U, 0xD8U, 0xDAU, 0xDDU, 0xE1U, 0xE4U, 0xE8U
};


}///////////////////////////////////////////////////////////////////////////////

template <>
struct prime< 1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <fixed_q Y_>
		XTAL_DEF_(return,inline,set)
		method()
		noexcept -> auto
		{
			return [&]<auto ...I> (bond::seek_t<I...>)
			XTAL_0FN_(to) (Y_{method(I)...})
				(bond::seek_s<fixed<Y_>::extent()>{});
		}

		template <int ...Ns_pow>
		XTAL_DEF_(return,inline,set)
		method(integral_q auto &&...is)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (0 == sizeof...(Ns_pow)) {return (one *...* method<1>(XTAL_REF_(is)));}
			XTAL_0IF (1 <= sizeof...(Ns_pow)) {return monomial_f<Ns_pow...>(method<1>(XTAL_REF_(is))...);}
		}
		template <int N_pow=1>
		XTAL_DEF_(return,inline,set)
		method(integral_variable_q auto i)
		noexcept -> auto
		{
			auto constexpr K_size = sizeof(_detail::U8_PRIME_HALF_SUBDEX) << 0;
			auto constexpr K_mask = K_size - one;

			assert(i < K_size);
			auto      o  = bond::math::bit_sign_f(i);
			i  &= ~o; o += 0 < i;
			i  &=  K_mask;
			i  += _detail::U8_PRIME_HALF_SUBDEX[i];
			i <<=  1;
			i  +=  o;
			return monomial_f<N_pow>(i);
		}

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(integral_constant_q auto i)
		noexcept -> auto
		{
			return constant_t<method<Ns...>(i())>{};
		}

	};
};
template <>
struct prime<-1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(integral_variable_q auto n)
		noexcept -> auto
		requires un_v<fixed_shaped_q<decltype(Ns)...>>
		{
			assert(0 <= n and n <= 719);
			auto constexpr *i0 = _std:: begin(_detail::U8_PRIME_HALF_SUBDEX);
			auto constexpr *iN = _std::   end(_detail::U8_PRIME_HALF_SUBDEX);
			auto const     *i_ = _std::lower_bound(i0, iN, n,
				[] (unsigned char const &x, extent_type const &y)
				XTAL_0FN_(to) (y > prime_f<1>(&x - i0)));

			return static_cast<XTAL_ALL_(n)>(i_ - i0);
		}

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(integral_constant_q auto i)
		noexcept -> auto
		requires un_v<fixed_shaped_q<decltype(Ns)...>>
		{
			return constant_t<method<Ns...>(i())>{};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <fixed_shaped_q auto M_nom>
struct prime<M_nom>
{
	XTAL_TYP_(set) M_typ = XTAL_ALL_(M_nom);
	XTAL_DEF_(set) M_ext = fixed_shaped<M_typ>::extent();
	XTAL_DEF_(set) M_num = prime_t<1>::template method<M_typ>();

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	protected:// OPERATE

		template <int N_ind, int N_dir=1>
		XTAL_DEF_(return,inline,set)
		factor_f(fixed_shaped_q auto &&x_)
		noexcept -> auto
		{
			using bond::math::bit_axis_f;
			auto const &x = get<N_ind>(x_);
			auto const &y = get<N_ind>(M_num);
			return exponential_f<>(_xtd::make_unsigned_f(bit_axis_f<N_dir>(x)), y);
		}

	public:// OPERATE

		template <integral_q auto N_ind>
		XTAL_DEF_(return,inline,set)
		edit(integral_q auto &n_product)
		noexcept -> auto
		{
			auto constexpr n_prime = prime_f<1>(N_ind);
			int i{};
			for (; 0 == n_product%n_prime; n_product /= n_prime) {++i;}
			return i;
		}

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method(fixed_shaped_q auto &&x_)
		noexcept -> auto
		{
			using          X_    =  XTAL_ALL_(x_);
			auto constexpr X_ext =  fixed_shaped<X_>::extent();
			static_assert(X_ext <= M_ext);

			auto const num = [&]<auto ...I> (bond::seek_t<I...>)
			XTAL_0FN_(to) (one *...* factor_f<I,  1>(x_))
				(bond::seek_s<X_ext>{});

			auto const nom = [&]<auto ...I> (bond::seek_t<I...>)
			XTAL_0FN_(to) (one *...* factor_f<I, -1>(x_))
				(bond::seek_s<X_ext>{});

			return bond::fit<X_>::ratio_f(num, nom);
		}
		/**/
		template <int ...Ns>
		XTAL_DEF_(return,inline,set)
		method(bond::seek_t<Ns...>)
		noexcept -> auto
		{
			return method(M_typ{Ns...});
		}
		/***/

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(real_q auto const u)
		noexcept -> auto
		{
			using U     = XTAL_ALL_(u);
			using U_fit = bond::fit<U>;
			auto  n = U_fit::delta_f(nearest_f<>(u*method<Ns...>(M_nom)));

			return [&]<auto ...I> (bond::seek_t<I...>)
			XTAL_0FN_(to) (M_typ{(edit<I>(n) - get<I>(M_nom))...})
				(bond::seek_s<M_ext>{});
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
