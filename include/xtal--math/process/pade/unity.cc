#pragma once
#include "./any.cc"
#include "./unity.hh"// testing...

#include "../dilated.hh"
#include "../dilate.hh"


XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_lim=0, int N_trim=0>
XTAL_DEF_(return,let)
unity_check_f(auto const &t)
noexcept -> bool
{
	int constexpr N_inf = -1;
	auto const u = unity_t<1>::template method_f<N_lim>(t);
	auto const v = unity_t<1>::template method_f<N_inf>(t);
	return check_f<N_trim>(u, v);
}


////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("unity")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using A_alpha = Eigen::ArrayXd ;// Eigen::Array<T_alpha,-1, 1>;
	using A_aphex = Eigen::ArrayXcd;// Eigen::Array<T_aphex,-1, 1>;

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	static_assert(same_q<T_alpha, decltype(U_phi{} (0))>);

	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	TRY_("flight")
	{
		T_alpha pitch = 0.1;
		T_alpha roll  = 0.1;
		T_alpha yaw   = 0.1;

		using T_foo = unity<1>;
		using W_foo = any_t<T_foo>;
		using Z_foo = typename W_foo::template dispatch<>;


		auto const w = unity_t<1>::template method_f<4>(pitch);
		auto const x = unity_t<1>::template method_f<4>(roll);
		auto const y = unity_t<1>::template method_f<4>(-yaw);

		auto const &[w_re, w_im] = destruct_f(w);
		auto const &[x_re, x_im] = destruct_f(x);
	//	auto const &[y_re, y_im] = destruct_f(y);

		T_aphex foo{w_re*x_re, w_im*x_im}; foo *= y;
		T_aphex bar{w_im*x_re, w_re*x_im}; bar *= y;

		atom::couple_t<T_alpha[4]> const o{foo.real(), bar.imag(), bar.real(), -foo.imag()};
		TRUE_(check_f<-2>(o[0], 0.73258330748146527));
		TRUE_(check_f<-2>(o[1], 0.10520194523965509));
		TRUE_(check_f<-2>(o[2], 0.66421893908191321));
		TRUE_(check_f<-2>(o[3], 0.10520194523965509));
	}
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

	// 1
		TRUE_(check_f<-40>(unity_t<1>::template method_f<-1>(x0), unity_t<1>::template method_f<1>(U_phi{x0, 0.})));
		TRUE_(check_f<-43>(unity_t<1>::template method_f<-1>(x1), unity_t<1>::template method_f<1>(U_phi{x1, 0.})));
		TRUE_(check_f<-43>(unity_t<1>::template method_f<-1>(x2), unity_t<1>::template method_f<1>(U_phi{x2, 0.})));
		TRUE_(check_f<-42>(unity_t<1>::template method_f<-1>(x3), unity_t<1>::template method_f<1>(U_phi{x3, 0.})));
		TRUE_(check_f<-42>(unity_t<1>::template method_f<-1>(x4), unity_t<1>::template method_f<1>(U_phi{x4, 0.})));
		TRUE_(check_f<-43>(unity_t<1>::template method_f<-1>(x5), unity_t<1>::template method_f<1>(U_phi{x5, 0.})));
		TRUE_(check_f<-43>(unity_t<1>::template method_f<-1>(x6), unity_t<1>::template method_f<1>(U_phi{x6, 0.})));
		TRUE_(check_f<-40>(unity_t<1>::template method_f<-1>(x7), unity_t<1>::template method_f<1>(U_phi{x7, 0.})));
	
	// 2
		TRUE_(check_f<-36>(unity_t<1>::template method_f<-1>(x0), unity_t<1>::template method_f<2>(U_phi{x0, 0.})));
		TRUE_(check_f<-37>(unity_t<1>::template method_f<-1>(x1), unity_t<1>::template method_f<2>(U_phi{x1, 0.})));
		TRUE_(check_f<-36>(unity_t<1>::template method_f<-1>(x2), unity_t<1>::template method_f<2>(U_phi{x2, 0.})));
		TRUE_(check_f<-36>(unity_t<1>::template method_f<-1>(x3), unity_t<1>::template method_f<2>(U_phi{x3, 0.})));
		TRUE_(check_f<-36>(unity_t<1>::template method_f<-1>(x4), unity_t<1>::template method_f<2>(U_phi{x4, 0.})));
		TRUE_(check_f<-36>(unity_t<1>::template method_f<-1>(x5), unity_t<1>::template method_f<2>(U_phi{x5, 0.})));
		TRUE_(check_f<-37>(unity_t<1>::template method_f<-1>(x6), unity_t<1>::template method_f<2>(U_phi{x6, 0.})));
		TRUE_(check_f<-36>(unity_t<1>::template method_f<-1>(x7), unity_t<1>::template method_f<2>(U_phi{x7, 0.})));
	
	//	3
		TRUE_(check_f<-27>(unity_t<1>::template method_f<-1>(x0), unity_t<1>::template method_f<3>(U_phi{x0, 0.})));
		TRUE_(check_f<-29>(unity_t<1>::template method_f<-1>(x1), unity_t<1>::template method_f<3>(U_phi{x1, 0.})));
		TRUE_(check_f<-27>(unity_t<1>::template method_f<-1>(x2), unity_t<1>::template method_f<3>(U_phi{x2, 0.})));
		TRUE_(check_f<-29>(unity_t<1>::template method_f<-1>(x3), unity_t<1>::template method_f<3>(U_phi{x3, 0.})));
		TRUE_(check_f<-29>(unity_t<1>::template method_f<-1>(x4), unity_t<1>::template method_f<3>(U_phi{x4, 0.})));
		TRUE_(check_f<-27>(unity_t<1>::template method_f<-1>(x5), unity_t<1>::template method_f<3>(U_phi{x5, 0.})));
		TRUE_(check_f<-29>(unity_t<1>::template method_f<-1>(x6), unity_t<1>::template method_f<3>(U_phi{x6, 0.})));
		TRUE_(check_f<-27>(unity_t<1>::template method_f<-1>(x7), unity_t<1>::template method_f<3>(U_phi{x7, 0.})));

	};
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
		A_aphex xs{{x0, x1, x2, x3, x4, x5, x6, x7}};

		auto y0 = unity_t<1>::template method_f<4>(x0);
		auto y1 = unity_t<1>::template method_f<4>(x1);
		auto y2 = unity_t<1>::template method_f<4>(x2);
		auto y3 = unity_t<1>::template method_f<4>(x3);
		auto y4 = unity_t<1>::template method_f<4>(x4);
		auto y5 = unity_t<1>::template method_f<4>(x5);
		auto y6 = unity_t<1>::template method_f<4>(x6);
		auto y7 = unity_t<1>::template method_f<4>(x7);
		auto ys = unity_t<1>::template method_f<4>(xs);

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

		auto y0 = unity_t<1>::template method_f<4>(x0);
		auto y1 = unity_t<1>::template method_f<4>(x1);
		auto y2 = unity_t<1>::template method_f<4>(x2);
		auto y3 = unity_t<1>::template method_f<4>(x3);
		auto y4 = unity_t<1>::template method_f<4>(x4);
		auto y5 = unity_t<1>::template method_f<4>(x5);
		auto y6 = unity_t<1>::template method_f<4>(x6);
		auto y7 = unity_t<1>::template method_f<4>(x7);
		auto ys = unity_t<1>::template method_f<4>(xs);

		TRUE_(y0.real(), ys.real(0));
		TRUE_(y1.real(), ys.real(1));
		TRUE_(y2.real(), ys.real(2));
		TRUE_(y3.real(), ys.real(3));
		TRUE_(y4.real(), ys.real(4));
		TRUE_(y5.real(), ys.real(5));
		TRUE_(y6.real(), ys.real(6));
		TRUE_(y7.real(), ys.real(7));

	};
	/***/
	TRY_("unity evaluation")
	{
		double const t0 = 1.125;

		TRUE_(unity_check_f<0,  2>(t0));
		TRUE_(unity_check_f<1,  6>(t0));
		TRUE_(unity_check_f<2, 14>(t0));
		TRUE_(unity_check_f<3, 20>(t0));
		TRUE_(unity_check_f<4, 28>(t0));
		TRUE_(unity_check_f<5, 37>(t0));
		TRUE_(unity_check_f<6, 47>(t0));
		TRUE_(unity_check_f<7, 49>(t0));

	}
	TRY_("dilution")
	{
		double const t1 = 1.11;
		double const t2 = 2.22;

		int constexpr N_lim = 4;

		auto z = unity_t<1>::template method_f<N_lim>(t1);

		TRUE_(check_f<-1>(z, process::lift_t<unity<1>, dilate<2>>::template method_f<N_lim>(t2)));
		TRUE_(check_f<-1>(z, process::lift_t<unity<1>           >::template method_f<N_lim>(t1)));

	}

	TRY_("unity<-1, -1> (native)")
	{
		T_aphex const x{0.123, 0.456};
		T_aphex const y = unity_t< 1>::template method_f<-1>(x);
		T_aphex const z = unity_t<-1>::template method_f<-1>(y);

		TRUE_(check_f<-19>(x, z));

	}
	TRY_("unity<-1,  5> (approx)")
	{
		T_aphex const x{0.123, 0.456};
		T_aphex const y = unity_t< 1>::template method_f< 5>(x);
		T_aphex const z = unity_t<-1>::template method_f< 5>(y);

		TRUE_(check_f<-19>(x, z));

	}
}
TAG_("unity trials")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using A_alpha = Eigen::ArrayXd ;// Eigen::Array<T_alpha,-1, 1>;
	using A_aphex = Eigen::ArrayXcd;// Eigen::Array<T_aphex,-1, 1>;

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	static_assert(same_q<T_alpha, decltype(U_phi{} (0))>);

	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	EST_("unity<+1; -1>\n   I^(4#)&\n   (*native floating-point*)")
	{
		T_alpha t = 0.618, _t = 0.414;
		T_aphex z {1, 0};

		for (T_sigma i = 0x100; ~--i;) {
			z *= unity_t<1>::template method_f<-1>(t);
			t += _t;
			t -= round(t);
		}
		return z;

	};
	EST_("unity<+1;  4>\n   I^(4#)&\n   (*approx floating-point*)")
	{
		T_alpha t = 0.618, _t = 0.414;
		T_aphex z {1, 0};

		for (T_sigma i = 0x100; ~--i;) {
			t += _t;
			t -= round(t);
			z *= unity_t<1>::template method_f<4>(t);
		}
		return z;

	};
	EST_("unity<+1;  4>\n   I^(4#)&\n   (*approx fixed-point*)")
	{
		U_phi t {0.618, 0.414};
		T_aphex z {1, 0};

		for (T_sigma i = 0x100; ~--i; ++t) {
			z *= unity_t<1>::template method_f<4>(t);
		}
		return z;

	};
	EST_("unity<-1, -1>\n   I^(4#)&\n   (*native floating-point*)")
	{
		T_aphex y{};
		for (int i = 0x20; ~--i;) {
			T_aphex const x{mt19937_f(), mt19937_f()};
			y += unity_t<-1>::template method_f<-1>(y);
		}
		return y;

	};
	EST_("unity<-1,  3>\n   I^(4#)&\n   (*approx floating-point*)")
	{
		T_aphex y{};
		for (int i = 0x20; ~--i;) {
			T_aphex const x{mt19937_f(), mt19937_f()};
			y += unity_t<-1>::template method_f< 3>(y);
		}
		return y;

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
