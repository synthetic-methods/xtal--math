#pragma once
#include "./any.cc"
#include "./tangy.hh"// testing...

#include "../dilated.hh"
#include "../dilate.hh"


XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("tangy")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using A_alpha = Eigen::Array<T_alpha,-1, 1>;
	using A_aphex = Eigen::Array<T_aphex,-1, 1>;

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	/**/
	TRY_("scalar evaluation")
	{
		TRUE_(check_f<51 - 20>(tangy_t< 2>::template method_f<-1>(T_aphex{1.0, 0.1}), tangy_t< 2>::template method_f< 6>(T_aphex{1.0, 0.1})));
		TRUE_(check_f<51 - 38>(tangy_t< 2>::template method_f<-1>(T_aphex{1.0, 0.1}), tangy_t< 2>::template method_f< 5>(T_aphex{1.0, 0.1})));
		TRUE_(check_f<51 - 32>(tangy_t< 2>::template method_f<-1>(T_aphex{1.0, 0.1}), tangy_t< 2>::template method_f< 4>(T_aphex{1.0, 0.1})));
		TRUE_(check_f<51 - 45>(tangy_t< 2>::template method_f<-1>(T_aphex{1.0, 0.1}), tangy_t< 2>::template method_f< 3>(T_aphex{1.0, 0.1})));
		TRUE_(check_f<51 - 41>(tangy_t< 2>::template method_f<-1>(T_aphex{1.0, 0.1}), tangy_t< 2>::template method_f< 2>(T_aphex{1.0, 0.1})));
	//	TRUE_(check_f<51 - 51>(tangy_t< 2>::template method_f<-1>(T_aphex{1.0, 0.1}), tangy_t< 2>::template method_f< 1>(T_aphex{1.0, 0.1})));

		TRUE_(check_f<-23>(tangy_t< 1>::template method_f<-1>(0.375), tangy_t< 1>::template method_f< 4>(0.375)));
		TRUE_(check_f<-38>(tangy_t< 1>::template method_f<-1>(0.375), tangy_t< 1>::template method_f< 3>(0.375)));
		TRUE_(check_f<-40>(tangy_t< 1>::template method_f<-1>(0.375), tangy_t< 1>::template method_f< 2>(0.375)));
		TRUE_(check_f<-50>(tangy_t< 1>::template method_f<-1>(0.375), tangy_t< 1>::template method_f< 1>(0.375)));
	//	TRUE_(check_f<-50>(tangy_t< 1>::template method_f<-1>(0.375), tangy_t< 1>::template method_f< 0>(0.375)));

		TRUE_(check_f<-23>(tangy_t< 1>::template method_f<-1>(0.125), tangy_t< 1>::template method_f< 4>(0.125)));
		TRUE_(check_f<-38>(tangy_t< 1>::template method_f<-1>(0.125), tangy_t< 1>::template method_f< 3>(0.125)));
		TRUE_(check_f<-40>(tangy_t< 1>::template method_f<-1>(0.125), tangy_t< 1>::template method_f< 2>(0.125)));
		TRUE_(check_f<-50>(tangy_t< 1>::template method_f<-1>(0.125), tangy_t< 1>::template method_f< 1>(0.125)));
	//	TRUE_(check_f<-50>(tangy_t< 1>::template method_f<-1>(0.125), tangy_t< 1>::template method_f< 0>(0.125)));

#if XTAL_ENV_(debug)
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template method_f<-1>(0.000), tangy_t<-1, 1>::template method_f< 5>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template method_f<-1>(0.000), tangy_t<-1, 1>::template method_f< 4>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template method_f<-1>(0.000), tangy_t<-1, 1>::template method_f< 3>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template method_f<-1>(0.000), tangy_t<-1, 1>::template method_f< 2>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template method_f<-1>(0.000), tangy_t<-1, 1>::template method_f< 1>(0.000)));
	//	TRUE_(check_f<-1>(tangy_t<-1, 1>::template method_f<-1>(0.000), tangy_t<-1, 1>::template method_f< 0>(0.000)));
#endif

	};
	/***/
	/**/
	TRY_("inverse evaluation")
	{
		T_aphex x;

		x = T_aphex{  0,  1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

		x = T_aphex{  0, -1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

		x = T_aphex{  1,  0};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

		x = T_aphex{ -1,  0};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

		x = T_aphex{  1,  1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

		x = T_aphex{  1, -1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

		x = T_aphex{ -1,  1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

		x = T_aphex{ -1, -1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));

	//	atan2(+<+)
		x = T_aphex{  2.414,  1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	//	atan2(+<-)
		x = T_aphex{ -2.414,  1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	//	atan2(-<+)
		x = T_aphex{  2.414, -1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	//	atan2(-<-)
		x = T_aphex{ -2.414, -1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	//	atan2(+>+)
		x = T_aphex{  1.618,  2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	//	atan2(+>-)
		x = T_aphex{ -1.618,  2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	//	atan2(->+)
		x = T_aphex{  1.618, -2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	//	atan2(->-)
		x = T_aphex{ -1.618, -2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_fit::patio_1, tangy_t<-1, 1>::template method_f< 7>(x.imag(), x.real())));

	};
	/***/

}
TAG_("tangy trials")
{
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using A_alpha = Eigen::Array<T_alpha,-1, 1>;
	using A_aphex = Eigen::Array<T_aphex,-1, 1>;

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	auto const mt19937_y = half*mt19937_f();

	EST_("tangy< 2; -1> (* Tanh *)\n~! native real")
	{
		return tangy_t< 2>::template method_f<-1>(mt19937_y);
	};
	EST_("tangy< 2;  3> (* Tanh *)\n~3 approx real")
	{
		return tangy_t< 2>::template method_f< 3>(mt19937_y);
	};
	EST_("tangy< 2;  2> (* Tanh *)\n~2 approx real")
	{
		return tangy_t< 2>::template method_f< 2>(mt19937_y);
	};
	EST_("tangy< 2;  1> (* Tanh *)\n~1 approx real")
	{
		return tangy_t< 2>::template method_f< 1>(mt19937_y);
	};
	EST_("tangy< 2;  0> (* Tanh *)\n~0 approx real")
	{
		return tangy_t< 2>::template method_f< 0>(mt19937_y);
	};


	EST_("tangy<-1; -1> (* ArcTan *)\n~! (reals)")
	{
		return tangy_t<-1>::template method_f<-1>(mt19937_y);
	};
	EST_("tangy<-1;  3> (* ArcTan *)\n~3 (reals)")
	{
		return tangy_t<-1>::template method_f< 3>(mt19937_y);
	};
	EST_("tangy<-1;  2> (* ArcTan *)\n~2 (reals)")
	{
		return tangy_t<-1>::template method_f< 2>(mt19937_y);
	};
	EST_("tangy<-1;  1> (* ArcTan *)\n~1 (reals)")
	{
		return tangy_t<-1>::template method_f< 1>(mt19937_y);
	};
	EST_("tangy<-1;  0> (* ArcTan *)\n~0 (reals)")
	{
		return tangy_t<-1>::template method_f< 0>(mt19937_y);
	};

}

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
