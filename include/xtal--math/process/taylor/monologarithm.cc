#pragma once
#include "./any.cc"

#include "../pade/unity.hh"



#include "./monologarithm.hh"
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("monologarithm")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;
	static constexpr U_alpha egg = 0.123456789;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	TRY_("evaluation")
	{
		TRUE_(check_f<  7>(monologarithm_t< 2>{}.template method<-1>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>{}.template method< 3>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>{}.template method< 2>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>{}.template method< 1>(egg), -log(1 - egg)));
		TRUE_(check_f<  7>(monologarithm_t< 2>{}.template method< 0>(egg), -log(1 - egg)));

		TRUE_(check_f<- 1>(monologarithm_t< 1>{}.template method<-1>(egg), +log(1 + egg)));
		TRUE_(check_f<- 2>(monologarithm_t< 1>{}.template method< 3>(egg), +log(1 + egg)));
		TRUE_(check_f<-13>(monologarithm_t< 1>{}.template method< 2>(egg), +log(1 + egg)));
		TRUE_(check_f< 18>(monologarithm_t< 1>{}.template method< 1>(egg), +log(1 + egg)));
		TRUE_(check_f< 10>(monologarithm_t< 1>{}.template method< 0>(egg), +log(1 + egg)));
		
		TRUE_(check_f<- 1>(monologarithm_t<-1>{}.template method<-1>(egg), exp(+egg) - 1));
		TRUE_(check_f<- 2>(monologarithm_t<-1>{}.template method< 3>(egg), exp(+egg) - 1));
		TRUE_(check_f<- 3>(monologarithm_t<-1>{}.template method< 2>(egg), exp(+egg) - 1));
		TRUE_(check_f< 24>(monologarithm_t<-1>{}.template method< 1>(egg), exp(+egg) - 1));
		TRUE_(check_f<  9>(monologarithm_t<-1>{}.template method< 0>(egg), exp(+egg) - 1));

		TRUE_(check_f<  7>(monologarithm_t<-2>{}.template method<-1>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>{}.template method< 3>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>{}.template method< 2>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>{}.template method< 1>(egg), 1 - exp(-egg)));
		TRUE_(check_f<  7>(monologarithm_t<-2>{}.template method< 0>(egg), 1 - exp(-egg)));

		TRUE_(check_f<24>(0.61803398874989490, monologarithm_t<-2,-1>{}.method(1.)));
		TRUE_(check_f<24>(1.61803398874989490, monologarithm_t<-1,-1>{}.method(1.)));
		TRUE_(check_f<24>(0.41421356237309515, monologarithm_t<-2,-1>{}.method(2.)));
		TRUE_(check_f<24>(2.41421356237309490, monologarithm_t<-1,-1>{}.method(2.)));
		TRUE_(check_f<24>(0.30277563773199456, monologarithm_t<-2,-1>{}.method(3.)));
		TRUE_(check_f<24>(3.30277563773199480, monologarithm_t<-1,-1>{}.method(3.)));

	}
	TRY_("mapping")
	{
		U_alpha zoom = U_fit::patio_f(1, 2);
		U_alpha s_abs = 0.88;
		U_alpha s_arg = 0.11;

		//\
		U_aphex w = pade::unity_t<1, dilate<-2>>{}.template method<4>(s_arg);
		//\
		U_aphex w = lift_t<pade::unity<1>, dilate<-2>>{}.template method<1>(s_arg);
		U_aphex w = link_t<pade::unity<1>, dilate<-2>>{}.template method<1>(s_arg);
		U_alpha u = one/s_abs;
		
		w *=      taylor::monologarithm_t<-1, 0>{}.template method<0>(u*zoom);
		w  = zoom/taylor::monologarithm_t<+1, 0>{}.template method<0>(w);
		w  = imagine_f<1>(w);

		TRUE_(check_f<16>(w, _std::complex{0.18009502457651236, 0.8570821020168073}));

	};
};
TAG_("monologarithm trials")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;
	static constexpr U_alpha egg = 0.123456789;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	auto mt19937_o = typename U_fit::MT19937{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	EST_("taylor::monologarithm_t<-1; -1>\n   Exp@# - 1&\n   (*native, floating-point*)")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = two*mt19937_f();
			//\
			w *= exp(x) - one;
			w *= monologarithm_t<-1>{}.template method<-1>(x);
		}
		return w;
	};
	EST_("taylor::monologarithm_t<-1;  2>\n   Exp@# - 1&\n   (*approx, floating-point*)")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = two*mt19937_f();
			w *= monologarithm_t<-1>{}.template method< 2>(x);
		}
		return w;
	};
	EST_("taylor::monologarithm_t<-1;  1>\n   Exp@# - 1&\n   (*approx, floating-point*)")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = two*mt19937_f();
			w *= monologarithm_t<-1>{}.template method< 1>(x);
		}
		return w;
	};
	EST_("taylor::monologarithm_t<-1;  0>\n   Exp@# - 1&\n   (*approx, floating-point*)")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = two*mt19937_f();
			w *= monologarithm_t<-1>{}.template method< 0>(x);
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
