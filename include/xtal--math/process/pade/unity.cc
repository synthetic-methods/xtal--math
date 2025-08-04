#pragma once
#include "./any.cc"
#include "./unity.hh"// testing...

#include "../dilated.hh"
#include "../dilate.hh"


XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("unity")
{
	using U_fit = bond::fit<>;

	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	using A_alpha = Eigen::ArrayXd ;// Eigen::Array<U_alpha,-1, 1>;
	using A_aphex = Eigen::ArrayXcd;// Eigen::Array<U_aphex,-1, 1>;

	static constexpr U_alpha two =  2;
	static constexpr U_alpha ten = 10;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	static_assert(same_q<U_alpha, decltype(U_phi{} (0))>);

	auto mt19937_o = typename U_fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	TRY_("unity<(+1)> // precision grading")
	{
		U_alpha x0{-0.4444444444444444};
		U_alpha x1{-0.3333333333333333};
		U_alpha x2{-0.2222222222222222};
		U_alpha x3{-0.1111111111111111};
		U_alpha x4{ 0.1111111111111111};
		U_alpha x5{ 0.2222222222222222};
		U_alpha x6{ 0.3333333333333333};
		U_alpha x7{ 0.4444444444444444};

	//	3 (- 1...- 1)
		TRUE_(check_f<- 1>(unity_t<(+1)>::template method_f<-1>(x0), unity_t<(+1)>::template method_f<3>(x0)));
		TRUE_(check_f<- 3>(unity_t<(+1)>::template method_f<-1>(x1), unity_t<(+1)>::template method_f<3>(x1)));
		TRUE_(check_f<- 1>(unity_t<(+1)>::template method_f<-1>(x2), unity_t<(+1)>::template method_f<3>(x2)));
		TRUE_(check_f<- 3>(unity_t<(+1)>::template method_f<-1>(x3), unity_t<(+1)>::template method_f<3>(x3)));
		TRUE_(check_f<- 3>(unity_t<(+1)>::template method_f<-1>(x4), unity_t<(+1)>::template method_f<3>(x4)));
		TRUE_(check_f<- 1>(unity_t<(+1)>::template method_f<-1>(x5), unity_t<(+1)>::template method_f<3>(x5)));
		TRUE_(check_f<- 3>(unity_t<(+1)>::template method_f<-1>(x6), unity_t<(+1)>::template method_f<3>(x6)));
		TRUE_(check_f<- 1>(unity_t<(+1)>::template method_f<-1>(x7), unity_t<(+1)>::template method_f<3>(x7)));

	//	2 (- 7...-12)
		TRUE_(check_f<- 7>(unity_t<(+1)>::template method_f<-1>(x0), unity_t<(+1)>::template method_f<2>(x0)));
		TRUE_(check_f<-12>(unity_t<(+1)>::template method_f<-1>(x1), unity_t<(+1)>::template method_f<2>(x1)));
		TRUE_(check_f<- 9>(unity_t<(+1)>::template method_f<-1>(x2), unity_t<(+1)>::template method_f<2>(x2)));
		TRUE_(check_f<-11>(unity_t<(+1)>::template method_f<-1>(x3), unity_t<(+1)>::template method_f<2>(x3)));
		TRUE_(check_f<-11>(unity_t<(+1)>::template method_f<-1>(x4), unity_t<(+1)>::template method_f<2>(x4)));
		TRUE_(check_f<- 9>(unity_t<(+1)>::template method_f<-1>(x5), unity_t<(+1)>::template method_f<2>(x5)));
		TRUE_(check_f<-12>(unity_t<(+1)>::template method_f<-1>(x6), unity_t<(+1)>::template method_f<2>(x6)));
		TRUE_(check_f<- 7>(unity_t<(+1)>::template method_f<-1>(x7), unity_t<(+1)>::template method_f<2>(x7)));

	//	1 (-19...-21)
		TRUE_(check_f<-19>(unity_t<(+1)>::template method_f<-1>(x0), unity_t<(+1)>::template method_f<1>(x0)));
		TRUE_(check_f<-21>(unity_t<(+1)>::template method_f<-1>(x1), unity_t<(+1)>::template method_f<1>(x1)));
		TRUE_(check_f<-16>(unity_t<(+1)>::template method_f<-1>(x2), unity_t<(+1)>::template method_f<1>(x2)));
		TRUE_(check_f<-21>(unity_t<(+1)>::template method_f<-1>(x3), unity_t<(+1)>::template method_f<1>(x3)));
		TRUE_(check_f<-21>(unity_t<(+1)>::template method_f<-1>(x4), unity_t<(+1)>::template method_f<1>(x4)));
		TRUE_(check_f<-16>(unity_t<(+1)>::template method_f<-1>(x5), unity_t<(+1)>::template method_f<1>(x5)));
		TRUE_(check_f<-21>(unity_t<(+1)>::template method_f<-1>(x6), unity_t<(+1)>::template method_f<1>(x6)));
		TRUE_(check_f<-19>(unity_t<(+1)>::template method_f<-1>(x7), unity_t<(+1)>::template method_f<1>(x7)));

	//	0 (-26...-29)
		TRUE_(check_f<-26>(unity_t<(+1)>::template method_f<-1>(x0), unity_t<(+1)>::template method_f<0>(x0)));
		TRUE_(check_f<-29>(unity_t<(+1)>::template method_f<-1>(x1), unity_t<(+1)>::template method_f<0>(x1)));
		TRUE_(check_f<-26>(unity_t<(+1)>::template method_f<-1>(x2), unity_t<(+1)>::template method_f<0>(x2)));
		TRUE_(check_f<-29>(unity_t<(+1)>::template method_f<-1>(x3), unity_t<(+1)>::template method_f<0>(x3)));
		TRUE_(check_f<-29>(unity_t<(+1)>::template method_f<-1>(x4), unity_t<(+1)>::template method_f<0>(x4)));
		TRUE_(check_f<-26>(unity_t<(+1)>::template method_f<-1>(x5), unity_t<(+1)>::template method_f<0>(x5)));
		TRUE_(check_f<-29>(unity_t<(+1)>::template method_f<-1>(x6), unity_t<(+1)>::template method_f<0>(x6)));
		TRUE_(check_f<-26>(unity_t<(+1)>::template method_f<-1>(x7), unity_t<(+1)>::template method_f<0>(x7)));

	};
	/**/
	TRY_("couple evaluation")
	{
		atom::couple_t<U_alpha[2]> x{0.123, 0.456};
		//\
		auto const v =        one/x;
		auto const v = U_alpha{1}/x;
		auto const y = unity_t<1>::template method_f<1>(x);

		TRUE_(check_f<8>(y.real(), atom::couple_t<U_alpha[2]>{ 0.715936483022, -0.962027671586}));
		TRUE_(check_f<8>(y.imag(), atom::couple_t<U_alpha[2]>{ 0.698165418993,  0.272951935517}));

	}
	/*/
	TRY_("vector evaluation")
	{
		U_aphex x0{ 0.0000000000000000, 0.0000000000000000};
		U_aphex x1{ 0.1111111111111111, 0.1111111111111111};
		U_aphex x2{ 0.2222222222222222, 0.2222222222222222};
		U_aphex x3{ 0.3333333333333333, 0.3333333333333333};
		U_aphex x4{ 0.4444444444444444, 0.4444444444444444};
		U_aphex x5{ 0.5555555555555555, 0.5555555555555555};
		U_aphex x6{ 0.6666666666666666, 0.6666666666666666};
		U_aphex x7{ 0.7777777777777777, 0.7777777777777777};
		A_aphex xs{{x0, x1, x2, x3, x4, x5, x6, x7}};

		auto y0 = unity_f<1>(x0);
		auto y1 = unity_f<1>(x1);
		auto y2 = unity_f<1>(x2);
		auto y3 = unity_f<1>(x3);
		auto y4 = unity_f<1>(x4);
		auto y5 = unity_f<1>(x5);
		auto y6 = unity_f<1>(x6);
		auto y7 = unity_f<1>(x7);
		auto ys = unity_f<1>(xs);

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
		U_aphex x0{ 0.0000000000000000, 0.0000000000000000};
		U_aphex x1{ 0.1111111111111111, 0.1111111111111111};
		U_aphex x2{ 0.2222222222222222, 0.2222222222222222};
		U_aphex x3{ 0.3333333333333333, 0.3333333333333333};
		U_aphex x4{ 0.4444444444444444, 0.4444444444444444};
		U_aphex x5{ 0.5555555555555555, 0.5555555555555555};
		U_aphex x6{ 0.6666666666666666, 0.6666666666666666};
		U_aphex x7{ 0.7777777777777777, 0.7777777777777777};
		U_alpha xs_co[2][8]
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

		auto y0 = unity_t<1>(x0);
		auto y1 = unity_t<1>(x1);
		auto y2 = unity_t<1>(x2);
		auto y3 = unity_t<1>(x3);
		auto y4 = unity_t<1>(x4);
		auto y5 = unity_t<1>(x5);
		auto y6 = unity_t<1>(x6);
		auto y7 = unity_t<1>(x7);
		auto ys = unity_t<1>(xs);

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

		TRUE_(check_f<-42>(unity_t<1, 1>::template method_f<(-1)>(t0), unity_t<1, 1>::template method_f<(-0)>(t0)));
		TRUE_(check_f<- 1>(unity_t<1, 1>::template method_f<(-1)>(t0), unity_t<1, 1>::template method_f<(-1)>(t0)));
		TRUE_(check_f<- 1>(unity_t<1, 1>::template method_f<(-1)>(t0), unity_t<1, 1>::template method_f<(-2)>(t0)));
		TRUE_(check_f<- 1>(unity_t<1, 1>::template method_f<(-1)>(t0), unity_t<1, 1>::template method_f<(-3)>(t0)));

	}
	TRY_("dilution")
	{
		double const t1 = 1.11;
		double const t2 = 2.22;

		int constexpr N_lim = 4;

		auto z = unity_t<1>::template method_f<N_lim>(t1);

		TRUE_(check_f<-1>(z, process::lift_t<unity<1, 1>, dilate<2>>::template method_f<N_lim>(t2)));
		TRUE_(check_f<-1>(z, process::lift_t<unity<1, 1>           >::template method_f<N_lim>(t1)));

	}

	TRY_("unity<-1, -1> (native)")
	{
		U_aphex const x{0.123, 0.456};
		U_aphex const y = unity_t< 1>::template method_f<-1>(x);
		U_aphex const z = unity_t<-1>::template method_f<-1>(y);

		TRUE_(check_f<-29>(x, z));

	}
	TRY_("unity<-1,  5> (approx)")
	{
		U_aphex const x{0.123, 0.456};
		U_aphex const y = unity_t< 1>::template method_f< 3>(x);
		U_aphex const z = unity_t<-1>::template method_f< 3>(y);

		TRUE_(check_f<-21>(x, z));

	}
	TRY_("unity<-1> corners")
	{
		U_aphex const x{0.123, 0.456};
		U_aphex const y = unity_t< 1>::template method_f< 3>(x);
		U_aphex const z = unity_t<-1>::template method_f< 3>(y);

		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{-2, 1}), U_aphex{ 0.42620819117478337, -0.12807499968169406}));
		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{-1, 2}), U_aphex{ 0.32379180882521663, -0.12807499968169406}));
		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{ 1, 2}), U_aphex{ 0.17620819117478337, -0.12807499968169406}));
		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{ 2, 1}), U_aphex{ 0.07379180882521663, -0.12807499968169406}));
		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{ 2,-1}), U_aphex{-0.07379180882521663, -0.12807499968169406}));
		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{ 1,-2}), U_aphex{-0.17620819117478337, -0.12807499968169406}));
		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{-1,-2}), U_aphex{-0.32379180882521663, -0.12807499968169406}));
		TRUE_(check_f<-24>(unity_t<-1>::template method_f< 3>(U_aphex{-2,-1}), U_aphex{-0.42620819117478337, -0.12807499968169406}));

	}
}
TAG_("unity trials")
{
	using U_fit = bond::fit<>;

	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	using A_alpha = Eigen::ArrayXd ;// Eigen::Array<U_alpha,-1, 1>;
	using A_aphex = Eigen::ArrayXcd;// Eigen::Array<U_aphex,-1, 1>;

	static constexpr U_alpha two =  2;
	static constexpr U_alpha ten = 10;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	static_assert(same_q<U_alpha, decltype(U_phi{} (0))>);

	auto mt19937_o = typename U_fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	EST_("unity<+1; -1>\n   I^(4#)&\n   (*native floating-point*)")
	{
		U_alpha t = 0.618, _t = 0.414;
		U_aphex z {1, 0};

		for (U_sigma i = 0x100; ~--i;) {
			z *= unity_t<1>::template method_f<-1>(t);
			t += _t;
			t -= round(t);
		}
		return z;

	};
	EST_("unity<+1; 1>\n   I^(4#)&\n   (*approx floating-point*)")
	{
		U_alpha t = 0.618, _t = 0.414;
		U_aphex z {1, 0};

		for (U_sigma i = 0x100; ~--i;) {
			t += _t;
			t -= round(t);
			z *= unity_t<1>::template method_f< 1>(t);
		}
		return z;

	};
	EST_("unity<+1; 1>\n   I^(4#)&\n   (*approx fixed-point*)")
	{
		U_phi t {0.618, 0.414};
		U_aphex z {1, 0};

		for (U_sigma i = 0x100; ~--i; ++t) {
			U_aphex const x{mt19937_f(), mt19937_f()};
			z *= unity_t<1>::template method_f< 1>(t);
		}
		return z;

	};
	EST_("unity<-1;-1>\n   I^(4#)&\n   (*native floating-point*)")
	{
		U_aphex y{};
		for (int i = 0x20; ~--i;) {
			U_aphex const x{mt19937_f(), mt19937_f()};
			y += unity_t<-1>::template method_f<-1>(x);
		}
		return y;

	};
	EST_("unity<-1; 1>\n   I^(4#)&\n   (*approx floating-point*)")
	{
		U_aphex y{};
		for (int i = 0x20; ~--i;) {
			U_aphex const x{mt19937_f(), mt19937_f()};
			y += unity_t<-1>::template method_f< 1>(x);
		}
		return y;

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
