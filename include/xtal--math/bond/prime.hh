#pragma once
#include "./any.hh"
#include "./bit.hh"





XTAL_ENV_(push)
namespace xtal::bond::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _detail
{///////////////////////////////////////////////////////////////////////////////

unsigned char constexpr primed[]
{	0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x02, 0x02
,	0x03, 0x05, 0x05, 0x07, 0x08, 0x08, 0x09, 0x0B
,	0x0D, 0x0D, 0x0F, 0x10, 0x10, 0x12, 0x13, 0x15
,	0x18, 0x19, 0x19, 0x1A, 0x1A, 0x1B, 0x21, 0x22
,	0x24, 0x24, 0x28, 0x28, 0x2A, 0x2C, 0x2D, 0x2F
,	0x31, 0x31, 0x35, 0x35, 0x36, 0x36, 0x3B, 0x40
,	0x41, 0x41, 0x42, 0x44, 0x44, 0x48, 0x4A, 0x4C
,	0x4E, 0x4E, 0x50, 0x51, 0x51, 0x55, 0x5B, 0x5C
,	0x5C, 0x5D, 0x63, 0x65, 0x69, 0x69, 0x6A, 0x6C
,	0x6F, 0x71, 0x73, 0x74, 0x76, 0x79, 0x7A, 0x7D
,	0x81, 0x81, 0x85, 0x85, 0x87, 0x88, 0x8A, 0x8D
,	0x8E, 0x8E, 0x8F, 0x94, 0x97, 0x98, 0x9B, 0x9C
,	0x9E, 0xA3, 0xA3, 0xAB, 0xAD, 0xB1, 0xB3, 0xB5
,	0xB5, 0xB7, 0xBB, 0xBD, 0xBF, 0xBF, 0xC1, 0xC3
,	0xC4, 0xC4, 0xC9, 0xCD, 0xCD, 0xCE, 0xD0, 0xD2
,	0xD2, 0xD7, 0xD8, 0xDA, 0xDD, 0xE1, 0xE4, 0xE8
};
XTAL_DEF_(return,inline,let)
primer_f(int n_index)
noexcept -> extent_type
{
	auto const i = n_index&0x7F; assert(i == n_index);
	return (static_cast<extent_type>(_detail::primed[i]) + i << 1)|!i + 1;
}
template <int N_index>
XTAL_DEF_(return,inline,let)
primer_f()
noexcept -> extent_type
{
	auto constexpr I = N_index&0x7F; static_assert(I == N_index);
	return primer_f(N_index);
}


}///////////////////////////////////////////////////////////////////////////////
/*!
\returns The prime number at the given index (valid for the first 128 primes).
*/
XTAL_DEF_(return,inline,let)
prime_f(int n_index)
noexcept -> extent_type
{
	auto const i = n_index&0x7F; assert(i == n_index);
	return (static_cast<extent_type>(_detail::primed[i]) + i << 1)|!i + 1;
}
XTAL_DEF_(return,inline,let)
prime_f(int n_index, int n_power)
noexcept -> extent_type
{
	n_power &= ~bit_sign_f(n_power);
	extent_type  m_value = 1;
	extent_type  n_value = prime_f(n_index);
	for (int x{n_power}; x; x >>= 1) {
		extent_type const w = -(1&x);
		m_value *= n_value&w|!w;
		n_value *= n_value;
	}
	return m_value;
}

template <int N_index>
XTAL_DEF_(return,inline,let)
prime_f(int n_power)
noexcept -> extent_type
{
	auto constexpr N = N_index&0x7F; static_assert(N == N_index);
	return prime_f(N, n_power);
}

template <int N_index, int N_power=1>
XTAL_DEF_(return,inline,let)
prime_f()
noexcept -> extent_type
{
	return prime_f<N_index>(N_power);
}

XTAL_DEF_(return,inline,let)
prime_f()
noexcept -> extent_type
{
	return 1;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!
\returns A prime root-of-unity-modulo-prime (valid for the first 50 primes).
*/
XTAL_DEF_(return,inline,let)
prime_recode_f(integral_q auto &n_product, int const n_index)
noexcept -> int
{
	auto const n_prime = prime_f(n_index);
	int i{};
	for (; 0 == n_product%n_prime; n_product /= n_prime) {++i;}
	return i;
}
XTAL_DEF_(return,inline,let)
prime_decode_f(integral_q auto n_product, int const n_index)
noexcept -> int
{
	return prime_recode_f(n_product, n_index);
}

template <int N_index, int N_scale=0>
XTAL_DEF_(return,inline,let)
prime_recode_f(integral_q auto &n_product)
noexcept -> int
{
	if constexpr (1 < N_scale) {
		return prime_recode_f(n_product, N_index) - prime_decode_f(N_scale, N_index);
	}
	else {
		return prime_recode_f(n_product, N_index);
	}
}
template <int N_index, int N_scale=0>
XTAL_DEF_(return,inline,let)
prime_decode_f(integral_q auto n_product)
noexcept -> int
{
	return prime_recode_f<N_index, N_scale>(n_product);
}



////////////////////////////////////////////////////////////////////////////////

template <unsigned N_shift, unsigned N_depth>
XTAL_DEF_(return,inline,let)
prime_encode_f()
noexcept -> auto
{
	return [&]<auto ...I> (bond::seek_t<I...>)
		XTAL_0FN_(to) (one *...* prime_f<I, N_shift>())
	(bond::seek_s<N_depth>{});
}

XTAL_DEF_(return,inline,let)
prime_encode_f(ordinal_q auto const ...i_primes)
noexcept -> auto
{
	auto constexpr I_pack = sizeof...(i_primes);
	auto const     i_pack = _std::tie(i_primes...);

	auto num = [&]<auto ...I> (bond::seek_t<I...>)
	XTAL_0FN_(to) (one *...* prime_f<I>( get<I>(i_pack)))
	(bond::seek_s<I_pack>{});

	auto nom = [&]<auto ...I> (bond::seek_t<I...>)
	XTAL_0FN_(to) (one *...* prime_f<I>(-get<I>(i_pack)))
	(bond::seek_s<I_pack>{});

	return num/nom;
}
template <int ...Ns>
XTAL_DEF_(return,inline,let)
prime_encode_f(bond::seek_t<Ns...>)
noexcept -> auto
{
	return prime_encode_f(Ns...);
}

template <unsigned N_scale>
XTAL_DEF_(return,inline,let)
prime_encode_f(ordinal_q auto const ...i_primes)
noexcept -> auto
requires (0 < sizeof...(i_primes))
{
	auto constexpr I_pack = sizeof...(i_primes);
	auto const     i_pack = _std::tie(i_primes...);

	return [&]<auto ...I> (bond::seek_t<I...>)
	XTAL_0FN_(to) (bond::fit<>::ratio_f(
		prime_encode_f((prime_decode_f<I>(N_scale) + get<I>(i_pack))...)
	,	N_scale
	))
	(bond::seek_s<I_pack>{});
}

template <unsigned N_shift, unsigned N_depth>
XTAL_DEF_(return,inline,let)
prime_encode_f(ordinal_q auto const ...i_primes)
noexcept -> auto
requires (0 < sizeof...(i_primes))
{
	auto constexpr N_scale = prime_encode_f<N_shift, N_depth>();
	return bond::fit<>::ratio_f(
		prime_encode_f((N_shift + i_primes)...)
	,	N_scale
	);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!
\returns A prime root-of-unity-modulo-prime (valid for the first 50 primes).
*/
XTAL_DEF_(return,inline,let)
prime_root_f(unsigned const n_prime, unsigned r_prime)
noexcept -> bool
{
	int o_count{0}, r_power{1};
	for (int i{0}; i < n_prime; ++i) {
		r_power *= r_prime;
		r_power %= n_prime;
		o_count += r_power == 1;
	}
	return o_count == 1 or n_prime == 2*r_prime;
}
XTAL_DEF_(return,inline,let)
prime_root_f(unsigned const n_prime)
noexcept -> unsigned
{
	assert(2 <= n_prime and n_prime <= 229 or n_prime == 257);
	auto p_prime = n_prime &  1;
	auto r_prime = p_prime +  1;
	auto m_prime = n_prime >> 1;
	m_prime &= (m_prime + 1)&(m_prime - 1);

	     if (3 <  n_prime    and not m_prime) {r_prime +=  1;}// => 5*7*17*31*127 * 257
	else if (0 == (  79*113*137*199)%n_prime) {r_prime +=  1;}// 28-bits
	else if (0 == ( 157*167*193*223)%n_prime) {r_prime +=  3;}// 31-bits
	else if (0 == (      97*151*229)%n_prime) {r_prime +=  5;}// 22-bits
	else if (0 == (23*71*73*103*109)%n_prime) {r_prime +=  9;}// 31-bits
	else if (0 == (41*43*47* 89*191)%n_prime) {r_prime += 17;}// 31-bits
	return r_prime;
}
template <unsigned N_prime>
XTAL_DEF_(return,inline,let)
prime_root_f()
noexcept -> unsigned
{
	static_assert(2 <= N_prime and N_prime <= 229 or N_prime == 257);
	return prime_root_f(N_prime);
}


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
