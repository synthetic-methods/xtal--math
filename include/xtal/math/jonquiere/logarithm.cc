#pragma once
#include "./any.cc"
#include "./logarithm.ii"// testing...






XTAL_ENV_(push)
namespace xtal::math::jonquiere::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("whatever")
{
	using re = bond::realized;

	using T_sigma = typename re::sigma_t;
	using T_delta = typename re::delta_t;
	using T_alpha = typename re::alpha_t;
	using T_aphex = typename re::aphex_t;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	using U_phi = group::cycle_t<T_alpha[2]>;

	auto mt19937_f = typename re::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("stuff")
	{
		TRUE_(check_f<7>(logarithm_t<+1, 0>::template function<-1>(0.123), -_std::log(1 - 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 0>::template function< 3>(0.123), -_std::log(1 - 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 0>::template function< 2>(0.123), -_std::log(1 - 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 0>::template function< 1>(0.123), -_std::log(1 - 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 0>::template function< 0>(0.123), -_std::log(1 - 0.123)));

		TRUE_(check_f<7>(logarithm_t<+1, 1>::template function<-1>(0.123), +_std::log(1 + 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 1>::template function< 3>(0.123), +_std::log(1 + 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 1>::template function< 2>(0.123), +_std::log(1 + 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 1>::template function< 1>(0.123), +_std::log(1 + 0.123)));
		TRUE_(check_f<7>(logarithm_t<+1, 1>::template function< 0>(0.123), +_std::log(1 + 0.123)));

		TRUE_(check_f<7>(logarithm_t<-1, 0>::template function<-1>(0.123), 1 - _std::exp(-0.123)));
		TRUE_(check_f<7>(logarithm_t<-1, 0>::template function< 3>(0.123), 1 - _std::exp(-0.123)));
		TRUE_(check_f<7>(logarithm_t<-1, 0>::template function< 2>(0.123), 1 - _std::exp(-0.123)));
		TRUE_(check_f<7>(logarithm_t<-1, 0>::template function< 1>(0.123), 1 - _std::exp(-0.123)));
		TRUE_(check_f<7>(logarithm_t<-1, 0>::template function< 0>(0.123), 1 - _std::exp(-0.123)));

		TRUE_(check_f<7>(logarithm_t<-1, 1>::template function<-1>(0.123), _std::exp(+0.123) - 1));
		TRUE_(check_f<7>(logarithm_t<-1, 1>::template function< 3>(0.123), _std::exp(+0.123) - 1));
		TRUE_(check_f<7>(logarithm_t<-1, 1>::template function< 2>(0.123), _std::exp(+0.123) - 1));
		TRUE_(check_f<7>(logarithm_t<-1, 1>::template function< 1>(0.123), _std::exp(+0.123) - 1));
		TRUE_(check_f<7>(logarithm_t<-1, 1>::template function< 0>(0.123), _std::exp(+0.123) - 1));

		TRUE_(check_f<24>(0.61803398874989490, logarithm_t<(-1), 0,-1>::template function<0>(1L)));
		TRUE_(check_f<24>(1.61803398874989490, logarithm_t<(-1), 1,-1>::template function<0>(1L)));
		TRUE_(check_f<24>(0.41421356237309515, logarithm_t<(-1), 0,-1>::template function<0>(2L)));
		TRUE_(check_f<24>(2.41421356237309490, logarithm_t<(-1), 1,-1>::template function<0>(2L)));
		TRUE_(check_f<24>(0.30277563773199456, logarithm_t<(-1), 0,-1>::template function<0>(3L)));
		TRUE_(check_f<24>(3.30277563773199480, logarithm_t<(-1), 1,-1>::template function<0>(3L)));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
