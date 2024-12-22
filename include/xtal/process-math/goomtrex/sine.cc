#pragma once
#include "./any.cc"
#include "./sine.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::goomtrex::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("sine")
{
	using _op = bond::operating;

	using     T_sigma = typename _op::sigma_type;
	using     T_delta = typename _op::delta_type;
	using     T_alpha = typename _op::alpha_type;
	using     T_aphex = typename _op::aphex_type;

	XTAL_SET_(T_alpha) one{1};
	XTAL_SET_(T_alpha) two{2};

	using U_phi = algebra::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	EST_("sine_t<2, -0>::function<~0>")
	{
		T_alpha o{};
		for (int i{0}; i < 0x100; ++i) {
			o += sinh(_op::mantissa_f(mt19937_f));
		}
		return o;
		
	};
	EST_("sine_t<2, -0>::function< 9>")
	{
		T_alpha o{};
		for (int i{0}; i < 0x100; ++i) {
			o += sine_t<2, -0>::template function<9>(_op::mantissa_f(mt19937_f));
		}
		return o;

	};
	EST_("sine_t<2, -0>::function< 3>")
	{
		T_alpha o{};
		for (int i{0}; i < 0x100; ++i) {
			o += sine_t<2, -0>::template function<9>(_op::mantissa_f(mt19937_f));
		}
		return o;

	};
	TRY_("sine_t<2, -0>")
	{
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<9>(one), 1.1742923611111113));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<8>(one), 1.1740493825137686));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<7>(one), 1.1736937830687832));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<6>(one), 1.1731430977472401));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<5>(one), 1.1675925925925927));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<4>(one), 1.1705016335204637));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<3>(one), 1.1666666666666667));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<2>(one), 1.1547005383792515));
		TRUE_(check_f<-1>(sine_t<2, -0>::template function<1>(one), 1.0000000000000000));

	};
	TRY_("sine_t<2, -1>")
	{
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<9>(one), 1.1742923611111113));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<8>(one), 1.1740493825137686));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<7>(one), 1.1736937830687832));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<6>(one), 1.1731430977472401));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<5>(one), 1.1675925925925927));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<4>(one), 1.1705016335204637));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<3>(one), 1.1666666666666667));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<2>(one), 1.1547005383792515));
		TRUE_(check_f<-1>(sine_t<2, -1>::template function<1>(one), 1.0000000000000000));

	};
	TRY_("sine_t<2, -2>")
	{
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<9>(one), 1.1742923611111113));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<8>(one), 1.1740493825137686));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<7>(one), 1.1736937830687832));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<6>(one), 1.1731430977472401));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<5>(one), 1.1675925925925927));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<4>(one), 1.1705016335204637));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<3>(one), 1.1666666666666667));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<2>(one), 1.1547005383792515));
		TRUE_(check_f<-1>(sine_t<2, -2>::template function<1>(one), 1.0000000000000000));

		TRUE_(true);

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
