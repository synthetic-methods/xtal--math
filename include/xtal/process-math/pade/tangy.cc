#pragma once
#include "./any.cc"
#include "./tangy.hh"// testing...

#include "../dilating.hh"
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

	using U_phi = algebra::d_::circular_t<T_alpha[2]>;

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

	};
	/***/
}

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
