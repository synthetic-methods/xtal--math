#pragma once
#include "./any.cc"
#include "./polylogarithm.hh"// testing...






XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_t;
	using T_delta = typename _op::delta_t;
	using T_alpha = typename _op::alpha_t;
	using T_aphex = typename _op::aphex_t;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	using U_phi = algebra::d_::circular_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		using _std::exp;
		using _std::log;

		TRUE_(check_f<7>(polylogarithm_t<+1>::template function<-1>(0.123), -log(1 - 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+1>::template function< 3>(0.123), -log(1 - 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+1>::template function< 2>(0.123), -log(1 - 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+1>::template function< 1>(0.123), -log(1 - 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+1>::template function< 0>(0.123), -log(1 - 0.123)));

		TRUE_(check_f<7>(polylogarithm_t<+2>::template function<-1>(0.123), +log(1 + 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+2>::template function< 3>(0.123), +log(1 + 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+2>::template function< 2>(0.123), +log(1 + 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+2>::template function< 1>(0.123), +log(1 + 0.123)));
		TRUE_(check_f<7>(polylogarithm_t<+2>::template function< 0>(0.123), +log(1 + 0.123)));

		TRUE_(check_f<7>(polylogarithm_t<-1>::template function<-1>(0.123), 1 - exp(-0.123)));
		TRUE_(check_f<7>(polylogarithm_t<-1>::template function< 3>(0.123), 1 - exp(-0.123)));
		TRUE_(check_f<7>(polylogarithm_t<-1>::template function< 2>(0.123), 1 - exp(-0.123)));
		TRUE_(check_f<7>(polylogarithm_t<-1>::template function< 1>(0.123), 1 - exp(-0.123)));
		TRUE_(check_f<7>(polylogarithm_t<-1>::template function< 0>(0.123), 1 - exp(-0.123)));

		TRUE_(check_f<7>(polylogarithm_t<-2>::template function<-1>(0.123), exp(+0.123) - 1));
		TRUE_(check_f<7>(polylogarithm_t<-2>::template function< 3>(0.123), exp(+0.123) - 1));
		TRUE_(check_f<7>(polylogarithm_t<-2>::template function< 2>(0.123), exp(+0.123) - 1));
		TRUE_(check_f<7>(polylogarithm_t<-2>::template function< 1>(0.123), exp(+0.123) - 1));
		TRUE_(check_f<7>(polylogarithm_t<-2>::template function< 0>(0.123), exp(+0.123) - 1));

		TRUE_(check_f<24>(0.61803398874989490, polylogarithm_f<-1,-1>(1L)));
		TRUE_(check_f<24>(1.61803398874989490, polylogarithm_f<-2,-1>(1L)));
		TRUE_(check_f<24>(0.41421356237309515, polylogarithm_f<-1,-1>(2L)));
		TRUE_(check_f<24>(2.41421356237309490, polylogarithm_f<-2,-1>(2L)));
		TRUE_(check_f<24>(0.30277563773199456, polylogarithm_f<-1,-1>(3L)));
		TRUE_(check_f<24>(3.30277563773199480, polylogarithm_f<-2,-1>(3L)));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
