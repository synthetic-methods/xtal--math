#pragma once
#include "./any.cc"
#include "./monologarithm.hh"// testing...

#include "../pade/unity.hh"



XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("monologarithm")
{
	using _op = bond::operating;

	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	using T_aphex = typename _op::aphex_type;
	static constexpr T_alpha egg =  0.123456789;
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	using U_phi = algebra::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _op::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("inversion")
	{
		TRUE_(check_f<-1>(one, monologarithm_t< 2,  1>::template function<0>(egg)*monologarithm_t< 2, -1>::template function<0>(egg)));
		TRUE_(check_f<-1>(one, monologarithm_t< 1,  1>::template function<0>(egg)*monologarithm_t< 1, -1>::template function<0>(egg)));
		TRUE_(check_f<-1>(one, monologarithm_t<-1,  1>::template function<0>(egg)*monologarithm_t<-1, -1>::template function<0>(egg)));
		TRUE_(check_f<-1>(one, monologarithm_t<-2,  1>::template function<0>(egg)*monologarithm_t<-2, -1>::template function<0>(egg)));

		TRUE_(check_f<19>(egg, monologarithm_t< 2,  1>::template function<0>(monologarithm_t<-2,  1>::template function<0>(egg))));
		TRUE_(check_f<19>(egg, monologarithm_t< 1,  1>::template function<0>(monologarithm_t<-1,  1>::template function<0>(egg))));
		TRUE_(check_f<19>(egg, monologarithm_t<-1,  1>::template function<0>(monologarithm_t< 1,  1>::template function<0>(egg))));
		TRUE_(check_f<19>(egg, monologarithm_t<-2,  1>::template function<0>(monologarithm_t< 2,  1>::template function<0>(egg))));

	}
	TRY_("evaluation")
	{
		TRUE_(check_f< 7>(monologarithm_t< 2>::template function<-1>(egg), -log(1 - egg)));
		TRUE_(check_f< 7>(monologarithm_t< 2>::template function< 3>(egg), -log(1 - egg)));
		TRUE_(check_f< 7>(monologarithm_t< 2>::template function< 2>(egg), -log(1 - egg)));
		TRUE_(check_f< 7>(monologarithm_t< 2>::template function< 1>(egg), -log(1 - egg)));
		TRUE_(check_f< 7>(monologarithm_t< 2>::template function< 0>(egg), -log(1 - egg)));

		TRUE_(check_f<- 1>(monologarithm_t< 1>::template function<-1>(egg), +log(1 + egg)));
		TRUE_(check_f<- 2>(monologarithm_t< 1>::template function< 3>(egg), +log(1 + egg)));
		TRUE_(check_f<-13>(monologarithm_t< 1>::template function< 2>(egg), +log(1 + egg)));
		TRUE_(check_f< 18>(monologarithm_t< 1>::template function< 1>(egg), +log(1 + egg)));
		TRUE_(check_f< 10>(monologarithm_t< 1>::template function< 0>(egg), +log(1 + egg)));
		
		TRUE_(check_f<-1>(monologarithm_t<-1>::template function<-1>(egg), exp(+egg) - 1));
		TRUE_(check_f<-2>(monologarithm_t<-1>::template function< 3>(egg), exp(+egg) - 1));
		TRUE_(check_f<-3>(monologarithm_t<-1>::template function< 2>(egg), exp(+egg) - 1));
		TRUE_(check_f<24>(monologarithm_t<-1>::template function< 1>(egg), exp(+egg) - 1));
		TRUE_(check_f< 9>(monologarithm_t<-1>::template function< 0>(egg), exp(+egg) - 1));

		TRUE_(check_f< 7>(monologarithm_t<-2>::template function<-1>(egg), 1 - exp(-egg)));
		TRUE_(check_f< 7>(monologarithm_t<-2>::template function< 3>(egg), 1 - exp(-egg)));
		TRUE_(check_f< 7>(monologarithm_t<-2>::template function< 2>(egg), 1 - exp(-egg)));
		TRUE_(check_f< 7>(monologarithm_t<-2>::template function< 1>(egg), 1 - exp(-egg)));
		TRUE_(check_f< 7>(monologarithm_t<-2>::template function< 0>(egg), 1 - exp(-egg)));

		TRUE_(check_f<24>(0.61803398874989490, monologarithm_f<-2, 1,-1>(1L)));
		TRUE_(check_f<24>(1.61803398874989490, monologarithm_f<-1, 1,-1>(1L)));
		TRUE_(check_f<24>(0.41421356237309515, monologarithm_f<-2, 1,-1>(2L)));
		TRUE_(check_f<24>(2.41421356237309490, monologarithm_f<-1, 1,-1>(2L)));
		TRUE_(check_f<24>(0.30277563773199456, monologarithm_f<-2, 1,-1>(3L)));
		TRUE_(check_f<24>(3.30277563773199480, monologarithm_f<-1, 1,-1>(3L)));

	}
	TRY_("mapping")
	{
		T_alpha zoom = _op::patio_f(1, 2);
		T_alpha s_abs = 0.88;
		T_alpha s_arg = 0.11;

		T_aphex w = pade::unity_t<1, dilate<1>>::template function<4>(s_arg);
		T_alpha u = _op::alpha_1/s_abs;
		
		w *= taylor::monologarithm_t<-1, +1>::template function<0>(u*zoom);
		w  = taylor::monologarithm_t<+1, -1>::template function<0>(w)*zoom;
		w  = imagine_f<1>(w);

		TRUE_(check_f<16>(w, _std::complex{0.18009502457651236, 0.8570821020168073}));

	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
