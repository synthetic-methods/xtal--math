#pragma once
#include "./any.cc"
#include "./symbol.hh"// testing...





XTAL_ENV_(push)
namespace xtal::algebra::math::jacobi::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("symbol")
{
	using _op = bond::operating;
	using T_delta = typename _op::delta_type;
	using T_sigma = typename _op::sigma_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;


	TRY_("11th characterization (integer)")
	{
		ordinal_type constexpr N = 11;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_delta[N]>;
		W w; w.characterize();

		TRUE_(w == W{0, 0, -1, 2, 3, 1,-4,-2,-3, 4,-5});

	}
	TRY_("11th subcharacterization (integer)")
	{
		ordinal_type constexpr N = 11;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_delta[K]>;
		W w; w.subcharacterize();

		TRUE_(w == W{0, -1, 2, 3, 1});

	}
	TRY_("11th characterization (complex)")
	{
		ordinal_type constexpr N = 11;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_aphex[N]>;
		W w; w.characterize();
		W m {_op::explo_f(w, K)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(_std::arg(XTAL_REF_(z))*K/_op::patio_1)>(w);

		_detail::apply_to<bond::computrim_f<16>>(m);
		_detail::apply_to<bond::computrim_f<16>>(w);

		TRUE_(m == W{0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1});
		TRUE_(w == W{0, 0, -1, 2, 3, 1,-4,-2,-3, 4,-5});

	}

	TRY_("7th characterization (integer)")
	{
		ordinal_type constexpr N = 7;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_delta[N]>;
		W w; w.characterize();

		TRUE_(w == W{ 0, 0, 2, 1,-2,-1,-3});

	}
	TRY_("7th subcharacterization (integer)")
	{
		ordinal_type constexpr N = 7;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_delta[K]>;
		W w; w.subcharacterize();

		TRUE_(w == W{ 0, 2, 1});

	}
	TRY_("7th characterization (complex)")
	{
		ordinal_type constexpr N = 7;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_aphex[N]>;
		W w; w.characterize();
		W m {_op::explo_f(w, K)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(_std::arg(XTAL_REF_(z))*K/_op::patio_1)>(w);

		_detail::apply_to<bond::computrim_f<16>>(m);
		_detail::apply_to<bond::computrim_f<16>>(w);

		TRUE_(m == W{0, 1, 1,-1, 1,-1,-1});
		TRUE_(w == W{0, 0, 2, 1,-2,-1,-3});

	}
	TRY_("7th subcharacterization (complex)")
	{
		ordinal_type constexpr N = 7;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_aphex[K]>;
		W w; w.subcharacterize();
		W m {_op::explo_f(w, K)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(_std::arg(XTAL_REF_(z))*K/_op::patio_1)>(w);

		_detail::apply_to<bond::computrim_f<16>>(m);
		_detail::apply_to<bond::computrim_f<16>>(w);

		TRUE_(w == W{0, 2, 1});

	}
	TRY_("7th characterization (real)")
	{
		ordinal_type constexpr N = 7;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_alpha[N]>;
		W w; w.characterize();

		_detail::apply_to<bond::computrim_f<16>>(w);

		TRUE_(w == W{0, 1, 1,-1, 1,-1,-1});

	}

	TRY_("5th characterization (complex)")
	{
		ordinal_type constexpr N = 5;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_aphex[N]>;
		W w; w.characterize();
		W m{_op::explo_f(w, K)};
		_detail::apply_to<[] (auto &&z) XTAL_0FN_(_std::arg(XTAL_REF_(z))*K/_op::patio_1)>(w);

		_detail::apply_to<bond::computrim_f<16>>(m);
		_detail::apply_to<bond::computrim_f<16>>(w);

		TRUE_(m == W{0, 1,-1,-1, 1});
		TRUE_(w == W{w[0], 0, 1,-1, w[M]});
	}
	/*/
	TRY_("5th characterization (real)")
	{
		ordinal_type constexpr N = 5;
		ordinal_type constexpr M = N  - 1;
		ordinal_type constexpr K = M >> 1;

		using W = symbol_t<T_alpha[N]>;
		W w; w.template characterize<2>();
		_detail::apply_to<bond::computrim_f<16>>(w);
	//	echo(w);
	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
