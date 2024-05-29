#pragma once
#include "./any.cc"
#include "./wnity.hh"// testing...

#include "../dilate.hh"





XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_lim=0, int N_trim=0>
XTAL_FN2 wnity__check_f(auto const &t)
XTAL_0EX
{
	using _std::get;

	int constexpr N_inf = -1;
	auto const u = wnity_t<1>::template function<N_lim>(t);
	auto const v = wnity_t<1>::template function<N_inf>(t);
	return true
		and bond::computrim_f<N_trim>(get<0>(u)) == bond::computrim_f<N_trim>(get<0>(v))
		and bond::computrim_f<N_trim>(get<1>(u)) == bond::computrim_f<N_trim>(get<1>(v))
		and true;
}


////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("wnity")
{
	using op = bond::operating;

	using T_sigma = typename op::sigma_t;
	using T_delta = typename op::delta_t;
	using T_alpha = typename op::alpha_t;
	using T_aphex = typename op::aphex_t;

	using A_alpha = Eigen::Array<T_alpha,-1, 1>;
	using A_aphex = Eigen::Array<T_aphex,-1, 1>;

	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	using U_phi = algebra::differential::circular_t<T_alpha[2]>;

	auto mt19937_f = typename op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	/*/
	TRY_("vector evaluation")
	{
		using _std::get;

		T_aphex x0{ 0.000000000000000, 0.000000000000000};
		T_aphex x1{ 0.111111111111111, 0.111111111111111};
		T_aphex x2{ 0.222222222222222, 0.222222222222222};
		T_aphex x3{ 0.333333333333333, 0.333333333333333};
		T_aphex x4{ 0.444444444444444, 0.444444444444444};
		T_aphex x5{ 0.555555555555555, 0.555555555555555};
		T_aphex x6{ 0.666666666666666, 0.666666666666666};
		T_aphex x7{ 0.777777777777777, 0.777777777777777};
		A_aphex xs{{x0, x1, x2, x3, x4, x5, x6, x7}};

		auto y0 = wnity_t<1>::template function<4>(x0);
		auto y1 = wnity_t<1>::template function<4>(x1);
		auto y2 = wnity_t<1>::template function<4>(x2);
		auto y3 = wnity_t<1>::template function<4>(x3);
		auto y4 = wnity_t<1>::template function<4>(x4);
		auto y5 = wnity_t<1>::template function<4>(x5);
		auto y6 = wnity_t<1>::template function<4>(x6);
		auto y7 = wnity_t<1>::template function<4>(x7);
		auto ys = wnity_t<1>::template function<4>(xs);
		auto [ys_0, ys_1] = ys;

		TRUE_(check_f<19>(y0.get<0>(), ys_0(0)));
		TRUE_(check_f<19>(y1.get<0>(), ys_0(1)));
		TRUE_(check_f<19>(y2.get<0>(), ys_0(2)));
		TRUE_(check_f<19>(y3.get<0>(), ys_0(3)));
		TRUE_(check_f<19>(y4.get<0>(), ys_0(4)));
		TRUE_(check_f<19>(y5.get<0>(), ys_0(5)));
		TRUE_(check_f<19>(y6.get<0>(), ys_0(6)));
		TRUE_(check_f<19>(y7.get<0>(), ys_0(7)));

	//	TODO: Compare performance of `processor` vs `Eigen`. \
	
	};
	/***/
	TRY_("evaluation (floating-point)")
	{
		double const t0 = 1.125;

		TRUE_(wnity__check_f<0,  2>(t0));
		TRUE_(wnity__check_f<1,  6>(t0));
		TRUE_(wnity__check_f<2, 14>(t0));
		TRUE_(wnity__check_f<3, 20>(t0));
		TRUE_(wnity__check_f<4, 28>(t0));
		TRUE_(wnity__check_f<5, 37>(t0));
		TRUE_(wnity__check_f<6, 47>(t0));
		TRUE_(wnity__check_f<7, 49>(t0));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
