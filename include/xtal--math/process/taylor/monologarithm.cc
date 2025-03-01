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
	using _fit = bond::fit<>;

	using T_sigma = typename _fit::sigma_type;
	using T_delta = typename _fit::delta_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;
	static constexpr T_alpha egg = 0.123456789;

	using U_phi = atom::math::phason_t<T_alpha[2]>;

	auto mt19937_f = typename _fit::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation")
	{
		TRUE_(check_f<  7>(monologarithm_t< 2>::template method_f<-1>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>::template method_f< 3>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>::template method_f< 2>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>::template method_f< 1>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>::template method_f< 0>(egg), -log(1 - egg)));

		TRUE_(check_f<- 1>(monologarithm_t< 1>::template method_f<-1>(egg), +log(1 + egg)));
		TRUE_(check_f<- 2>(monologarithm_t< 1>::template method_f< 3>(egg), +log(1 + egg)));
		TRUE_(check_f<-13>(monologarithm_t< 1>::template method_f< 2>(egg), +log(1 + egg)));
		TRUE_(check_f< 18>(monologarithm_t< 1>::template method_f< 1>(egg), +log(1 + egg)));
		TRUE_(check_f< 10>(monologarithm_t< 1>::template method_f< 0>(egg), +log(1 + egg)));
		
		TRUE_(check_f<- 1>(monologarithm_t<-1>::template method_f<-1>(egg), exp(+egg) - 1));
		TRUE_(check_f<- 2>(monologarithm_t<-1>::template method_f< 3>(egg), exp(+egg) - 1));
		TRUE_(check_f<- 3>(monologarithm_t<-1>::template method_f< 2>(egg), exp(+egg) - 1));
		TRUE_(check_f< 24>(monologarithm_t<-1>::template method_f< 1>(egg), exp(+egg) - 1));
		TRUE_(check_f<  9>(monologarithm_t<-1>::template method_f< 0>(egg), exp(+egg) - 1));

		TRUE_(check_f<  7>(monologarithm_t<-2>::template method_f<-1>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>::template method_f< 3>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>::template method_f< 2>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>::template method_f< 1>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>::template method_f< 0>(egg), 1 - exp(-egg)));

		TRUE_(check_f<24>(0.61803398874989490, monologarithm_t<-2,-1>::method_f(1.)));
		TRUE_(check_f<24>(1.61803398874989490, monologarithm_t<-1,-1>::method_f(1.)));
		TRUE_(check_f<24>(0.41421356237309515, monologarithm_t<-2,-1>::method_f(2.)));
		TRUE_(check_f<24>(2.41421356237309490, monologarithm_t<-1,-1>::method_f(2.)));
		TRUE_(check_f<24>(0.30277563773199456, monologarithm_t<-2,-1>::method_f(3.)));
		TRUE_(check_f<24>(3.30277563773199480, monologarithm_t<-1,-1>::method_f(3.)));

	}
	TRY_("mapping")
	{
		T_alpha zoom = _fit::patio_f(1, 2);
		T_alpha s_abs = 0.88;
		T_alpha s_arg = 0.11;

		//\
		T_aphex w = pade::unity_t<1, dilate<2>>::template method_f<4>(s_arg);
		T_aphex w = lift_t<pade::unity<1>, dilate<2>>::template method_f<4>(s_arg);
		T_alpha u = one/s_abs;
		
		w *=      taylor::monologarithm_t<-1, 0>::template method_f<0>(u*zoom);
		w  = zoom/taylor::monologarithm_t<+1, 0>::template method_f<0>(w);
		w  = imagine_f<1>(w);

		TRUE_(check_f<16>(w, _std::complex{0.18009502457651236, 0.8570821020168073}));

	};
	EST_("real std::exp(...) - 1")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) - one;
			w *= exp(x) - one;
		}
		return w;
	};
	EST_("real antimonologarithm... <M_iso=-2, M_lim=-1>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) - one;
			w *= monologarithm_t<-2>::template method_f<~0>(x);
		}
		return w;
	};
	EST_("real antimonologarithm... <M_iso=-2, M_lim= 0>")
	{
		T_alpha w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto x = _fit::mantissa_f(mt19937_f) - one;
			w *= monologarithm_t<-2>::template method_f< 0>(x);
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
