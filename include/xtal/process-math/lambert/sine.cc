#pragma once
#include "./any.cc"
#include "./sine.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::lambert::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("sine< 1>")
{
	using _op = bond::operate<size_type>;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("task")
	{
		double constexpr N_2pi = 6.283185307179586476925286766559005768L;
		int constexpr N_lim = 3;
		int constexpr N_bas = 0;// e
		
		double const t = 0.5, t_2pi = N_2pi*t;

		TRUE_(check_f<-1>(0.5210953054937474, sine_t< 2, -0>::template function<~0>(0.5)));
		TRUE_(check_f<-1>(0.5207289109917556, sine_t< 2, -0>::template function< 3>(0.5)));
		TRUE_(check_f<-1>(0.5196377064333049, sine_t< 2, -0>::template function< 2>(0.5)));
		TRUE_(check_f<-1>(0.5153882032022076, sine_t< 2, -0>::template function< 1>(0.5)));

		TRUE_(check_f<-50>(std:: sinh(2.)/(1.), sine_t< 2, -0>::template function<3>(2.)));
		TRUE_(check_f<-50>(std:: sinh(2.)/(2.), sine_t< 2, -1>::template function<3>(2.)));

	}
	EST_("real std::sinh")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += _std::sinh(x);
		}
		return w;
	};
	EST_("real area sine... <N_lim=3>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += sine_t< 2, 0>::template function<3>(x);
		}
		return w;
	};
	EST_("real area sine... <N_lim=2>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += sine_t< 2, 0>::template function<2>(x);
		}
		return w;
	};
	EST_("real area sine... <N_lim=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += sine_t< 2, 0>::template function<1>(x);
		}
		return w;
	};
}
/***/
/**/
TAG_("sine< 1>")
{
	using _op = bond::operate<size_type>;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("task")
	{
		double constexpr N_2pi = 6.283185307179586476925286766559005768L;
		int constexpr N_lim = 3;
		int constexpr N_bas = 0;// e
		
		double const t = 0.5, t_2pi = N_2pi*t;

		TRUE_(check_f<-1>(0.4812118250596035, sine_t<-2, -0>::template function<~0>(0.5)));
		TRUE_(check_f<-4>(0.4815020643585153, sine_t<-2, -0>::template function< 3>(0.5)));
		TRUE_(check_f<-1>(0.4823734124990402, sine_t<-2, -0>::template function< 2>(0.5)));
		TRUE_(check_f<-1>(0.4858682717566457, sine_t<-2, -0>::template function< 1>(0.5)));

		TRUE_(check_f<-46>(std::asinh(2.)/(1.), sine_t<-2, -0>::template function<3>(2.)));
		TRUE_(check_f<-46>(std::asinh(2.)/(2.), sine_t<-2, -1>::template function<3>(2.)));

	}
	EST_("real std::asinh")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += _std::asinh(x);
		}
		return w;
	};
	EST_("real area sine... <N_lim=3>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += sine_t<-2, 0>::template function<3>(x);
		}
		return w;
	};
	EST_("real area sine... <N_lim=2>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += sine_t<-2, 0>::template function<2>(x);
		}
		return w;
	};
	EST_("real area sine... <N_lim=1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _op::mantissa_f(mt19937_f) + one;
			w += sine_t<-2, 0>::template function<1>(x);
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
