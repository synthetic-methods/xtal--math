#pragma once
#include "./any.cc"
#include "./symbol.hh"// testing...

#include "../../process/all.hh"



XTAL_ENV_(push)
namespace xtal::atom::math::dirichlet::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("symbol")
{
	using _fit = bond::fit<>;
	using T_delta = typename _fit::delta_type;
	using T_sigma = typename _fit::sigma_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;


	TRY_("11th characterization (integer)")
	{
		int constexpr N = 11;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_delta[N]>;
		W w; w.characterize();

		TRUE_(w == W{0, 0, -1, 2, 3, 1,-4,-2,-3, 4,-5});

	}
	TRY_("11th subcharacterization (integer)")
	{
		int constexpr N = 11;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_delta[K]>;
		W w; w.subcharacterize();

		TRUE_(w == W{0, -1, 2, 3, 1});

	}
	TRY_("11th characterization (complex)")
	{
		int constexpr N = 11;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_aphex[N]>;
		W w; w.characterize();
		W m {process::math::power_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(w);

		TRUE_(m == W{0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1});
		TRUE_(w == W{0, 0, -1, 2, 3, 1,-4,-2,-3, 4,-5});

	}

	TRY_("7th characterization (integer)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_delta[N]>;
		W w; w.characterize();

		TRUE_(w == W{ 0, 0, 2, 1,-2,-1,-3});

	}
	TRY_("7th subcharacterization (integer)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_delta[K]>;
		W w; w.subcharacterize();

		TRUE_(w == W{ 0, 2, 1});

	}
	TRY_("7th characterization (complex)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_aphex[N]>;
		W w; w.characterize();
		W m {process::math::power_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(w);

		TRUE_(m == W{0, 1, 1,-1, 1,-1,-1});
		TRUE_(w == W{0, 0, 2, 1,-2,-1,-3});

	}
	TRY_("7th subcharacterization (complex)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_aphex[K]>;
		W w; w.subcharacterize();
		W m {process::math::power_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(w);

		TRUE_(w == W{0, 2, 1});

	}
	TRY_("7th characterization (real)")
	{
		int constexpr N = 7;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_alpha[N]>;
		W w; w.characterize();

		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(w);

		TRUE_(w == W{0, 1, 1,-1, 1,-1,-1});

	}

	TRY_("5th characterization (complex)")
	{
		int constexpr N = 5;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_aphex[N]>;
		W w; w.characterize();
		W m{process::math::power_f<K>(w)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(to) (arg(XTAL_REF_(z))*K/_fit::patio_1)>(w);

		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(m);
		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(w);

		TRUE_(m == W{0, 1,-1,-1, 1});
		TRUE_(w == W{w[0], 0, 1,-1, w[M]});
	}
	/*/
	TRY_("5th characterization (real)")
	{
		int constexpr N = 5;
		int constexpr M = N  - 1;
		int constexpr K = M >> 1;

		using W = symbol_t<T_alpha[N]>;
		W w; w.template characterize<2>();
		_detail::apply_to<[] XTAL_1FN_(function) (_fit::template trim_f<16>)>(w);
	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
