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
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation")
	{
		TRUE_(check_f<-35>(one, cbrt(-0.500)*root_f<-3>(-0.500)));
		TRUE_(check_f<-34>(one, cbrt(-0.250)*root_f<-3>(-0.250)));
		TRUE_(check_f<-33>(one, cbrt(-0.125)*root_f<-3>(-0.125)));

	//	TRUE_(check_f<-1>(1.0f, cbrt(0.500f)*root_f<-3>(0.500f)));
	//	TRUE_(check_f<-1>(1.0f, cbrt(0.250f)*root_f<-3>(0.250f)));
	//	TRUE_(check_f<-7>(1.0f, cbrt(0.125f)*root_f<-3>(0.125f)));
	//	NOTE: GCC aborts on `TRUE_` (though the `check_f` itself succeeds)...

		TRUE_(check_f<-1>(root_f< 4>(0.5),     sqrt(sqrt(0.5))));
		TRUE_(check_f<-1>(root_f<-4>(0.5), one/sqrt(sqrt(0.5))));

		TRUE_(check_f<-1>(root_f< 8>(0.5),     sqrt(sqrt(sqrt(0.5)))));
		TRUE_(check_f<-1>(root_f<-8>(0.5), one/sqrt(sqrt(sqrt(0.5)))));

		TRUE_(check_f<22>(root_f<-2>(pow(T_aphex {2, 3}, -2.0)), T_aphex {2, 3}));
		TRUE_(check_f<22>(pow(T_aphex {2, 3}, 0.5), root_f< 2>(T_aphex {2, 3})));

		TRUE_(check_f<-1>(one/sqrt(2.2345268795805384), root_f<-2>(_std::complex{2.2345268795805384,0.0}).real()));

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
			auto x = _op::mantissa_f(mt19937_f) + one;
			w *= _std::cbrt(x);
		}
		return w;
	};
	EST_("real root_t<-3>::<1>")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w *= root_t<-3>::template function<1>(x);
		}
		return w;
	};
	EST_("real root_t<-3>::<2>")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w *= root_t<-3>::template function<2>(x);
		}
		return w;
	};
	EST_("real root_t<-3>::<3>")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w *= root_t<-3>::template function<3>(x);
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
