#pragma once
#include "./any.cc"

#include "../../process/all.hh"



#include "./symbol.hh"
XTAL_ENV_(push)
namespace xtal::atom::math::dirichlet::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("symbol")
{
	using U_fit   = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	TRY_("11th characterization (integer)")
	{
		int constexpr N = 11;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_delta[N]>;
		W w; w.characterize();

		TRUE_(w == W{0, 0, 1,-2, 2, 4,-1,-3, 3,-4,-5});

	}
	TRY_("11th subcharacterization (integer)")
	{
		int constexpr N = 11;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_delta[K]>;

		W w; w.subcharacterize();

		TRUE_(w == W{0, -1, 2, 3, 1});

	}
	TRY_("11th characterization (complex)")
	{
		int constexpr N = 11;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_aphex[N]>;
		W w; w.characterize();
		W m {process::math::monomial_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/U_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(w);

		TRUE_(m == W{0, 1,-1, 1, 1, 1,-1,-1,-1, 1,-1});
		TRUE_(w == W{0, 0, 1,-2, 2, 4,-1,-3, 3,-4,-5});

	}

	TRY_("7th characterization (integer)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_delta[N]>;
		W w; w.characterize();

		TRUE_(w == W{ 0, 0, 2, 1,-2,-1,-3});

	}
	TRY_("7th subcharacterization (integer)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_delta[K]>;
		W w; w.subcharacterize();

		TRUE_(w == W{ 0, 2, 1});

	}
	TRY_("7th characterization (complex)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_aphex[N]>;
		W w; w.characterize();
		W m {process::math::monomial_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/U_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(w);

		TRUE_(m == W{0, 1, 1,-1, 1,-1,-1});
		TRUE_(w == W{0, 0, 2, 1,-2,-1,-3});

	}
	TRY_("7th subcharacterization (complex)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_aphex[K]>;
		W w; w.subcharacterize();
		W m {process::math::monomial_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/U_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(w);

		TRUE_(w == W{0, 2, 1});

	}
	TRY_("7th characterization (real)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_alpha[N]>;
		W w; w.characterize();

		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(w);

		TRUE_(w == W{0, 1, 1,-1, 1,-1,-1});

	}

	TRY_("5th characterization (complex)")
	{
		int constexpr N = 5;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_aphex[N]>;
		W w; w.characterize();
		W m{process::math::monomial_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/U_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(w);

		TRUE_(m == W{0, 1,-1,-1, 1});
	//	TRUE_(w == W{w[0], 0, 1,-1, w[M]});
	}
	/*/
	TRY_("5th characterization (real)")
	{
		int constexpr N = 5;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<U_alpha[N]>;
		W w; w.template characterize<2>();
		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(w);
	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
