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
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;

	using A_alpha = Eigen::Array<T_alpha,-1, 1>;
	using A_aphex = Eigen::Array<T_aphex,-1, 1>;

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = arrange::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	/**/
	TRY_("scalar evaluation")
	{
		TRUE_(check_f<-23>(tangy_t< 1>::template function<-1>(0.375), tangy_t< 1>::template function< 4>(0.375)));
		TRUE_(check_f<-37>(tangy_t< 1>::template function<-1>(0.375), tangy_t< 1>::template function< 3>(0.375)));
		TRUE_(check_f<-40>(tangy_t< 1>::template function<-1>(0.375), tangy_t< 1>::template function< 2>(0.375)));
		TRUE_(check_f<-50>(tangy_t< 1>::template function<-1>(0.375), tangy_t< 1>::template function< 1>(0.375)));
		TRUE_(check_f<-50>(tangy_t< 1>::template function<-1>(0.375), tangy_t< 1>::template function< 0>(0.375)));

		TRUE_(check_f<-23>(tangy_t< 1>::template function<-1>(0.125), tangy_t< 1>::template function< 4>(0.125)));
		TRUE_(check_f<-37>(tangy_t< 1>::template function<-1>(0.125), tangy_t< 1>::template function< 3>(0.125)));
		TRUE_(check_f<-40>(tangy_t< 1>::template function<-1>(0.125), tangy_t< 1>::template function< 2>(0.125)));
		TRUE_(check_f<-50>(tangy_t< 1>::template function<-1>(0.125), tangy_t< 1>::template function< 1>(0.125)));
		TRUE_(check_f<-50>(tangy_t< 1>::template function<-1>(0.125), tangy_t< 1>::template function< 0>(0.125)));

#if XTAL_ENV_(debug)
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template function<-1>(0.000), tangy_t<-1, 1>::template function< 5>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template function<-1>(0.000), tangy_t<-1, 1>::template function< 4>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template function<-1>(0.000), tangy_t<-1, 1>::template function< 3>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template function<-1>(0.000), tangy_t<-1, 1>::template function< 2>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template function<-1>(0.000), tangy_t<-1, 1>::template function< 1>(0.000)));
		TRUE_(check_f<-1>(tangy_t<-1, 1>::template function<-1>(0.000), tangy_t<-1, 1>::template function< 0>(0.000)));
#endif
	};
	/***/
	/**/
	TRY_("inverse evaluation")
	{
		T_aphex x;

		x = T_aphex{  0,  1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

		x = T_aphex{  0, -1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

		x = T_aphex{  1,  0};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

		x = T_aphex{ -1,  0};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

		x = T_aphex{  1,  1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

		x = T_aphex{  1, -1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

		x = T_aphex{ -1,  1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

		x = T_aphex{ -1, -1};
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-1>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));

	//	atan2(+<+)
		x = T_aphex{  2.414,  1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	//	atan2(+<-)
		x = T_aphex{ -2.414,  1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	//	atan2(-<+)
		x = T_aphex{  2.414, -1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	//	atan2(-<-)
		x = T_aphex{ -2.414, -1.618};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	//	atan2(+>+)
		x = T_aphex{  1.618,  2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	//	atan2(+>-)
		x = T_aphex{ -1.618,  2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	//	atan2(->+)
		x = T_aphex{  1.618, -2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	//	atan2(->-)
		x = T_aphex{ -1.618, -2.414};
		TRUE_(check_f<-40>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 0>(x.imag(), x.real())));
		TRUE_(check_f<-35>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 1>(x.imag(), x.real())));
		TRUE_(check_f<-29>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 2>(x.imag(), x.real())));
		TRUE_(check_f<-26>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 3>(x.imag(), x.real())));
		TRUE_(check_f<-21>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 4>(x.imag(), x.real())));
		TRUE_(check_f<-16>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 5>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 6>(x.imag(), x.real())));
		TRUE_(check_f<-13>(atan2(x.imag(), x.real())/_op::patio_1, tangy_t<-1, 1>::template function< 7>(x.imag(), x.real())));

	};
	/***/

	EST_("tangy inversion<-1> (native)")
	{
		T_aphex y{};
		for (int i = 0x20; ~--i;) {
			T_aphex const x{_op::mantissa_f(mt19937_f), _op::mantissa_f(mt19937_f)};
			y += tangy_t<-1, 1>::template function<-1>(x.imag(), x.real());
		}
		return y;

	};
	EST_("tangy inversion< 4> (approximation)")
	{
		T_aphex y{};
		for (int i = 0x20; ~--i;) {
			T_aphex const x{_op::mantissa_f(mt19937_f), _op::mantissa_f(mt19937_f)};
			y += tangy_t<-1, 1>::template function< 4>(x.imag(), x.real());
		}
		return y;

	};
//	EST_("tangy inversion< 2> (approximation)")
//	{
//		T_aphex y{};
//		for (int i = 0x20; ~--i;) {
//			T_aphex const x{_op::mantissa_f(mt19937_f), _op::mantissa_f(mt19937_f)};
//			y += tangy_t<-1, 1>::template function< 2>(x.imag(), x.real());
//		}
//		return y;
//
//	};
//	EST_("tangy inversion< 1> (approximation)")
//	{
//		T_aphex y{};
//		for (int i = 0x20; ~--i;) {
//			T_aphex const x{_op::mantissa_f(mt19937_f), _op::mantissa_f(mt19937_f)};
//			y += tangy_t<-1, 1>::template function< 1>(x.imag(), x.real());
//		}
//		return y;
//
//	};
//	EST_("tangy inversion< 0> (approximation)")
//	{
//		T_aphex y{};
//		for (int i = 0x20; ~--i;) {
//			T_aphex const x{_op::mantissa_f(mt19937_f), _op::mantissa_f(mt19937_f)};
//			y += tangy_t<-1, 1>::template function< 0>(x.imag(), x.real());
//		}
//		return y;
//
//	};

}

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
