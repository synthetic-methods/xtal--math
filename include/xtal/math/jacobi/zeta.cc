#pragma once
#include "./any.cc"
#include "./zeta.hh"// testing...






XTAL_ENV_(push)
namespace xtal::math::jacobi::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("whatever")
{
	using op = bond::operating;

	using U_sigma = typename op::sigma_t;
	using U_delta = typename op::delta_t;
	using U_alpha = typename op::alpha_t;
	using U_aphex = typename op::aphex_t;

	using W_alpha = Eigen::Array<U_alpha, -1, 1>;
	using W_aphex = Eigen::Array<U_aphex, -1, 1>;

	//\
//	static_assert(multiplicative_group_p<2, W_aphex>);
//	static_assert(equality_p<2, W_aphex>);
	static_assert(complex_field_q<W_aphex>);

	auto mt19937_f = typename op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("scalar evaluation")
	{
		using V_phi = algebra::differential::circular_t<U_alpha[2]>;
		using U_phi = _std::complex<V_phi>;

		//\
		U_phi   constexpr t = { 0.22, 0.11};
		U_aphex constexpr t = { 0.22, 0.11};
		U_aphex constexpr s = { 0.22, 1.11};
		U_aphex constexpr r = {-0.150796, 0.106422};

		TRUE_(bond::computrim_f<16>(r) == bond::computrim_f<16>(zeta_f(t, s)));

	};
	/*/
	TRY_("vector evaluation")
	{
		W_aphex t{{ 0.22, 0.11}};
		W_aphex s{{ 0.22, 1.11}};
		W_aphex r{{-0.150796, 0.106422}};

		auto _r = zeta_f(t, s);

	//	TRUE_(bond::computrim_f<16>(r) == bond::computrim_f<16>(zeta_f(t, s)));

	};
	/***/
}

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
