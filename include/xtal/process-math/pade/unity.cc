#pragma once
#include "./any.cc"
#include "./unity.hh"// testing...

#include "../dilate.hh"
#include "../dilute.hh"




XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_lim=0, int N_trim=0>
XTAL_FN2 unity__check_f(auto const &t)
XTAL_0EX
{
	int constexpr N_inf = -1;
	auto const u = unity_t<1>::template function<N_lim>(t);
	auto const v = unity_t<1>::template function<N_inf>(t);
	return bond::computrim_f<N_trim>(u) == bond::computrim_f<N_trim>(v);
}


////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("unity")
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

	using U_phi = algebra::d_::circular_t<T_alpha[2]>;

	auto mt19937_f = typename op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	/**/
	TRY_("scalar evaluation")
	{
		T_alpha x0{-0.4444444444444444};
		T_alpha x1{-0.3333333333333333};
		T_alpha x2{-0.2222222222222222};
		T_alpha x3{-0.1111111111111111};
		T_alpha x4{ 0.1111111111111111};
		T_alpha x5{ 0.2222222222222222};
		T_alpha x6{ 0.3333333333333333};
		T_alpha x7{ 0.4444444444444444};

		TRUE_(check_f<18>(unity_t<1>::template function<4>(x0), unity_t<1>::template function<4>(U_phi{x0, 0.})));
		TRUE_(check_f<18>(unity_t<1>::template function<4>(x1), unity_t<1>::template function<4>(U_phi{x1, 0.})));
		TRUE_(check_f<18>(unity_t<1>::template function<4>(x2), unity_t<1>::template function<4>(U_phi{x2, 0.})));
		TRUE_(check_f<18>(unity_t<1>::template function<4>(x3), unity_t<1>::template function<4>(U_phi{x3, 0.})));
		TRUE_(check_f<18>(unity_t<1>::template function<4>(x4), unity_t<1>::template function<4>(U_phi{x4, 0.})));
		TRUE_(check_f<18>(unity_t<1>::template function<4>(x5), unity_t<1>::template function<4>(U_phi{x5, 0.})));
		TRUE_(check_f<18>(unity_t<1>::template function<4>(x6), unity_t<1>::template function<4>(U_phi{x6, 0.})));
		TRUE_(check_f<18>(unity_t<1>::template function<4>(x7), unity_t<1>::template function<4>(U_phi{x7, 0.})));

	};
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
		A_aphex xs{{x0, x1, x2, x3, x4, x5, x6, x7}};

	//	unity_t<1> _y1{}; _y1 <<= V_unity_limit(3); echo(_y1(x1));
		auto y0 = unity_t<1>::template function<4>(x0);
		auto y1 = unity_t<1>::template function<4>(x1);
		auto y2 = unity_t<1>::template function<4>(x2);
		auto y3 = unity_t<1>::template function<4>(x3);
		auto y4 = unity_t<1>::template function<4>(x4);
		auto y5 = unity_t<1>::template function<4>(x5);
		auto y6 = unity_t<1>::template function<4>(x6);
		auto y7 = unity_t<1>::template function<4>(x7);
		auto ys = unity_t<1>::template function<4>(xs);

		TRUE_(check_f<19>(y0, ys(0)));
		TRUE_(check_f<19>(y1, ys(1)));
		TRUE_(check_f<19>(y2, ys(2)));
		TRUE_(check_f<19>(y3, ys(3)));
		TRUE_(check_f<19>(y4, ys(4)));
		TRUE_(check_f<19>(y5, ys(5)));
		TRUE_(check_f<19>(y6, ys(6)));
		TRUE_(check_f<19>(y7, ys(7)));

	};
	/***/
	TRY_("evaluation (floating-point)")
	{
		double const t0 = 1.125;

		TRUE_(unity__check_f<0,  2>(t0));
		TRUE_(unity__check_f<1,  6>(t0));
		TRUE_(unity__check_f<2, 14>(t0));
		TRUE_(unity__check_f<3, 20>(t0));
		TRUE_(unity__check_f<4, 28>(t0));
		TRUE_(unity__check_f<5, 37>(t0));
		TRUE_(unity__check_f<6, 47>(t0));
		TRUE_(unity__check_f<7, 49>(t0));

		EST_("evaluation <N_lim=-1>")
		{
			using _std::round;

			T_alpha t = 0.618, _t = 0.414;
			T_aphex z {1, 0};

			for (T_sigma i = 0x100; ~--i;) {
				z *= unity_t<1>::template function<-1>(t);
				t += _t;
				t -= round(t);
			}
			return z;

		};
		EST_("bench <N_lim=4>")
		{
			using _std::round;

			T_alpha t = 0.618, _t = 0.414;
			T_aphex z {1, 0};

			for (T_sigma i = 0x100; ~--i;) {
				t += _t;
				t -= round(t);
				z *= unity_t<1>::template function<4>(t);
			}
			return z;

		};
	}
	TRY_("evaluation (fixed-point)")
	{
		//\
		auto const [t0, t1] = algebra::d_::circular_t<T_alpha[2]> {1.125, 0.0};
		T_alpha t0{1.125};

		TRUE_(unity__check_f<0,  2>(t0));
		TRUE_(unity__check_f<1,  6>(t0));
		TRUE_(unity__check_f<2, 14>(t0));
		TRUE_(unity__check_f<3, 20>(t0));
		TRUE_(unity__check_f<4, 28>(t0));
		TRUE_(unity__check_f<5, 37>(t0));
		TRUE_(unity__check_f<6, 47>(t0));
		TRUE_(unity__check_f<7, 49>(t0));

		EST_("bench <N_lim=4>")
		{
			U_phi t {0.618, 0.414};
			T_aphex z {1, 0};

			for (T_sigma i = 0x100; ~--i; ++t) {
				z *= unity_t<1>::template function<4>(t);
			}
			return z;

		};
	}
	TRY_("dilution")
	{
		double const t1 = 1.11;
		double const t2 = 2.22;

		int constexpr N_lim = 4;

		auto z = unity_t<1>::template function<N_lim>(t1);

		TRUE_(check_f<-1>(z, unity_t<1           >::template function<N_lim>(t1)));
		TRUE_(check_f<-1>(z, unity_t<1, dilute<1>>::template function<N_lim>(t2)));
		
		TRUE_(check_f<-1>(z, process::chain_t<unity<1>, dilute<1>>::template function<N_lim>(t2)));
		TRUE_(check_f<-1>(z, process::chain_t<unity<1>           >::template function<N_lim>(t1)));


	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
