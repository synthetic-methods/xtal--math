#pragma once
#include "./any.cc"
#include "./root.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("root")
{
	using _op = bond::operating;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	XTAL_LET M_1i_root2 = root_f< 2>(T_aphex{1, 1});
	XTAL_LET M_i1_root2 = root_f<-2>(T_aphex{1, 1});

	XTAL_LET N_zero_root2 = root_f< 2>(_op::alpha_0);
	XTAL_LET N_half_root2 = root_f<-2>(_op::diplo_1);
	XTAL_LET N_half_root3 = root_f<-3>(_op::diplo_1);
	XTAL_LET N_half_root4 = root_f<-4>(_op::diplo_1);
	XTAL_LET N_half_root5 = root_f<-5>(_op::diplo_1);

	TRY_("evaluation")
	{
		TRUE_(check_f<-1>(root_f<-2>(T_aphex{1, 1}), one/root_f< 2>(T_aphex{1, 1})));

		TRUE_(check_f<-1>(N_half_root2, 0.7071067811865475244008443621048490L));
		TRUE_(check_f<-1>(N_half_root3, 0.7937005259840997373758528196361541L));
		TRUE_(check_f<-1>(N_half_root4, 0.8408964152537145430311254762332149L));
		TRUE_(check_f<-1>(N_half_root5, 0.8705505632961241391362700174797461L));
	//	echo(check_f(pow(0.5, _op::ratio_f(1,  3)), root_t< 3>::template function<2>(0.5)));
	//	echo(check_f(pow(0.5, _op::ratio_f(1, -3)), root_t<-3>::template function<2>(0.5)));

	//	echo(check_f(pow(0.5, _op::ratio_f(1,  3)), root_t< 3>::template function<3>(0.5)));
	//	echo(check_f(pow(0.5, _op::ratio_f(1, -3)), root_t<-3>::template function<3>(0.5)));

		TRUE_(check_f<-27>(pow(0.5, _op::ratio_f(1,  3)), root_f< 3>(0.5)));
		TRUE_(check_f<-28>(pow(0.5, _op::ratio_f(1, -3)), root_f<-3>(0.5)));

		TRUE_(check_f<-30>(pow(0.5, _op::ratio_f(1,  5)), root_f< 5>(0.5)));
		TRUE_(check_f<-27>(pow(0.5, _op::ratio_f(1, -5)), root_f<-5>(0.5)));

		TRUE_(check_f<-30>(pow(0.5, _op::ratio_f(1,  7)), root_f< 7>(0.5)));
		TRUE_(check_f<-28>(pow(0.5, _op::ratio_f(1, -7)), root_f<-7>(0.5)));

		TRUE_(check_f<-35>(1.0, cbrt(-0.500)*root_f<-3>(-0.500)));
		TRUE_(check_f<-34>(1.0, cbrt(-0.250)*root_f<-3>(-0.250)));
		TRUE_(check_f<-33>(1.0, cbrt(-0.125)*root_f<-3>(-0.125)));

	//	TRUE_(check_f<-1>(1.0f, cbrt(0.500f)*root_f<-3>(0.500f)));
	//	TRUE_(check_f<-1>(1.0f, cbrt(0.250f)*root_f<-3>(0.250f)));
	//	TRUE_(check_f<-7>(1.0f, cbrt(0.125f)*root_f<-3>(0.125f)));
	//	NOTE: GCC aborts on `TRUE_` (though the `check_f` itself succeeds)...

		TRUE_(check_f<-1>(root_f< 4>(0.5),     sqrt(sqrt(0.5))));
		TRUE_(check_f<-1>(root_f<-4>(0.5), 1.0/sqrt(sqrt(0.5))));

		TRUE_(check_f<-1>(root_f< 8>(0.5),     sqrt(sqrt(sqrt(0.5)))));
		TRUE_(check_f<-1>(root_f<-8>(0.5), 1.0/sqrt(sqrt(sqrt(0.5)))));

		TRUE_(check_f<22>(root_f<-2>(pow(T_aphex {2, 3}, -2.0)), T_aphex {2, 3}));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 0.5), root_f< 2>(T_aphex {2, 3})));

		TRUE_(check_f<-1>(1.0/sqrt(2.2345268795805384), root_f<-2>(_std::complex{2.2345268795805384,0.0}).real()));

	}
	TRY_("punctured evaluation")
	{
		TRUE_(0.0 == root_f<-1, 1>(0.0) - root_f<-1, 1>(0.0));
		TRUE_(0.0 == root_f<-2, 1>(0.0) - root_f<-2, 1>(0.0));

	}

	EST_("real 1/std::cbrt")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= one/cbrt(_op::mantissa_f(mt19937_f) + one);
		}
		return w;
	};
	EST_("real root_t<-3>::<3>")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= root_t<-3>::template function<3>(_op::mantissa_f(mt19937_f) + one);
		}
		return w;
	};
	EST_("real root_t<-3>::<2>")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= root_t<-3>::template function<2>(_op::mantissa_f(mt19937_f) + one);
		}
		return w;
	};
	EST_("real root_t<-3>::<1>")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= root_t<-3>::template function<1>(_op::mantissa_f(mt19937_f) + one);
		}
		return w;
	};

}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
