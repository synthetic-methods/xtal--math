#pragma once
#include "./any.cc"





#include "./proot.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test::XTAL_NUM
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
using namespace xtal::process::math::_test;

XTAL_DEF_(return,inline,let)
check_proot_f(unsigned const n_prime, unsigned r_prime)
noexcept -> bool
{
	unsigned o_count{0}, r_power{1};
	for (unsigned i{0}; i < n_prime; ++i) {
		r_power *= r_prime;
		r_power %= n_prime;
		o_count += r_power == 1U;
	}
	return o_count == 1U or n_prime == (r_prime << 1);
}


////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("proot")
{
	using _fit = bond::fit<>;
	using U_sigma = typename _fit::sigma_type;
	using U_delta = typename _fit::delta_type;
	using U_alpha = typename _fit::alpha_type;
	using U_aphex = typename _fit::aphex_type;

	TRY_("proot_f<1>(...)")
	{
		#pragma unroll
		for (int i{1}; i < 0x80; ++i) {
			TRUE_(check_proot_f(prime_f<1>(i), proot_f<1>(i)));
		}
	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
