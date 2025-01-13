#pragma once
#include "./any.cc"
#include "./wnity.hh"// testing...

#include "../dilated.hh"



XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_lim=0, int N_trim=0>
XTAL_DEF_(return)
XTAL_LET wnity_check_f(auto const &t)
noexcept -> bool
{
	int constexpr N_inf = -1;
	auto const u = wnity_t<1>::template function<N_lim>(t);
	auto const v = wnity_t<1>::template function<N_inf>(t);
	return true
		and check_f<N_trim>(get<0>(u), get<0>(v))
		and check_f<N_trim>(get<1>(u), get<1>(v))
		and true;
}


////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("wnity")
{
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;

	using A_alpha = Eigen::Array<T_alpha,-1, 1>;
	using A_aphex = Eigen::Array<T_aphex,-1, 1>;

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = arrange::math::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	/*/
	TRY_("vector evaluation")
	{
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

		TRUE_(check_f<19>(get<0>(y0), ys_0(0)));
		TRUE_(check_f<19>(get<0>(y1), ys_0(1)));
		TRUE_(check_f<19>(get<0>(y2), ys_0(2)));
		TRUE_(check_f<19>(get<0>(y3), ys_0(3)));
		TRUE_(check_f<19>(get<0>(y4), ys_0(4)));
		TRUE_(check_f<19>(get<0>(y5), ys_0(5)));
		TRUE_(check_f<19>(get<0>(y6), ys_0(6)));
		TRUE_(check_f<19>(get<0>(y7), ys_0(7)));

	//	TODO: Compare performance of `processor` vs `Eigen`. \
	
	};
	/***/
	/*/
	TRY_("vector evaluation")
	{
		T_aphex x0{ 0.0000000000000000, 0.0000000000000000};
		T_aphex x1{ 0.1111111111111111, 0.1111111111111111};
		T_aphex x2{ 0.2222222222222222, 0.2222222222222222};
		T_aphex x3{ 0.3333333333333333, 0.3333333333333333};
		T_aphex x4{ 0.4444444444444444, 0.4444444444444444};
		T_aphex x5{ 0.5555555555555555, 0.5555555555555555};
		T_aphex x6{ 0.6666666666666666, 0.6666666666666666};
		T_aphex x7{ 0.7777777777777777, 0.7777777777777777};
		T_alpha xs_co[2][8]
		{	{	0.0000000000000000,
				0.1111111111111111,
				0.2222222222222222,
				0.3333333333333333,
				0.4444444444444444,
				0.5555555555555555,
				0.6666666666666666,
				0.7777777777777777
			},
			{	0.0000000000000000,
				0.1111111111111111,
				0.2222222222222222,
				0.3333333333333333,
				0.4444444444444444,
				0.5555555555555555,
				0.6666666666666666,
				0.7777777777777777
			}
		};
		Eigen::Map<A_alpha> xs_re(xs_co[0], 8);
		Eigen::Map<A_alpha> xs_im(xs_co[1], 8);

		_std::complex<A_alpha> xs{xs_re, xs_im};

		auto y0 = wnity_t<1>::template function<4>(x0);
		auto y1 = wnity_t<1>::template function<4>(x1);
		auto y2 = wnity_t<1>::template function<4>(x2);
		auto y3 = wnity_t<1>::template function<4>(x3);
		auto y4 = wnity_t<1>::template function<4>(x4);
		auto y5 = wnity_t<1>::template function<4>(x5);
		auto y6 = wnity_t<1>::template function<4>(x6);
		auto y7 = wnity_t<1>::template function<4>(x7);
		auto ys = wnity_t<1>::template function<4>(xs);

		TRUE_(get<0>(y0).real(), get<0>(ys).real(0));
		TRUE_(get<0>(y1).real(), get<0>(ys).real(1));
		TRUE_(get<0>(y2).real(), get<0>(ys).real(2));
		TRUE_(get<0>(y3).real(), get<0>(ys).real(3));
		TRUE_(get<0>(y4).real(), get<0>(ys).real(4));
		TRUE_(get<0>(y5).real(), get<0>(ys).real(5));
		TRUE_(get<0>(y6).real(), get<0>(ys).real(6));
		TRUE_(get<0>(y7).real(), get<0>(ys).real(7));
	
	};
	/***/
	TRY_("evaluation (floating-point)")
	{
		double const t0 = 1.125;

		TRUE_(wnity_check_f<0,  2>(t0));
		TRUE_(wnity_check_f<1,  6>(t0));
		TRUE_(wnity_check_f<2, 14>(t0));
		TRUE_(wnity_check_f<3, 20>(t0));
		TRUE_(wnity_check_f<4, 28>(t0));
		TRUE_(wnity_check_f<5, 37>(t0));
		TRUE_(wnity_check_f<6, 47>(t0));
		TRUE_(wnity_check_f<7, 49>(t0));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
