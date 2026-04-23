#pragma once
#include "./any.cc"

#include "./tangent.hh"



#include "./logarithm.hh"
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("logarithm")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;
	static constexpr U_alpha egg =  1.23456789;
	static constexpr U_alpha ten = 10;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	/**/
	TRY_("evaluation")
	{
		TRUE_(check_f<- 8>(logarithm_t< 1, 1>{}.template method< 3>(egg), log(egg)));
		TRUE_(check_f<- 2>(logarithm_t<-1, 1>{}.template method< 3>(egg), exp(egg)));

		TRUE_(check_f<- 1>(logarithm_t< 1   >{}.template method<-1>(egg), log(egg)));
		TRUE_(check_f<- 1>(logarithm_t< 1   >{}.template method< 3>(egg), log(egg)));
		TRUE_(check_f<- 5>(logarithm_t< 1   >{}.template method< 2>(egg), log(egg)));
		TRUE_(check_f<-19>(logarithm_t< 1   >{}.template method< 1>(egg), log(egg)));
		TRUE_(check_f<-40>(logarithm_t< 1   >{}.template method< 0>(egg), log(egg)));

		TRUE_(check_f<- 1>(logarithm_t<-1   >{}.template method<-1>(egg), exp(egg)));
		TRUE_(check_f<- 1>(logarithm_t<-1   >{}.template method< 3>(egg), exp(egg)));
		TRUE_(check_f<-13>(logarithm_t<-1   >{}.template method< 2>(egg), exp(egg)));
		TRUE_(check_f<-29>(logarithm_t<-1   >{}.template method< 1>(egg), exp(egg)));
	//	UNTRUE_(check_f(logarithm_t<-1>{}.template method< 0>(egg), exp(egg)));

		TRUE_(check_f<-48>(logarithm_t< 1, 1>{}.template method<~0>(U_aphex{0.3, 0.8}), logarithm_t< 1, 1>{}.template method< 3>(U_aphex{0.3, 0.8})));
		TRUE_(check_f<-48>(logarithm_t< 1, 1>{}.template method<~0>(U_aphex{0.3, 0.8}), logarithm_t< 1, 1>{}.template method< 0>(U_aphex{0.3, 0.8})));

	}
	/***/
}
TAG_("logarithm trials")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	auto mt19937_o = typename U_fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	EST_("logarithm_t<+1,  1; -1>\n   Log@#&\n   (*complex, native*)")
	{
		U_aphex w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			auto y = mt19937_f() + 1.0;
			w += log(U_aphex{x, y});
		//	w += logarithm_t< 1, 1>{}.template method<-1>(U_aphex{x, y});
		}
		return w;
	};
	EST_("logarithm_t<+1,  1;  2>\n   Log@#&\n   (*complex, approx*)")
	{
		U_aphex w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			auto y = mt19937_f() + 1.0;
			w += logarithm_t< 1, 1>{}.template method< 2>(U_aphex{x, y});
		}
		return w;
	};
	EST_("logarithm_t<+1,  0;  2>\n   Log@#&\n   (*complex, approx*)")
	{
		U_aphex w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			auto y = mt19937_f() + 1.0;
			w += logarithm_t< 1, 0>{}.template method< 2>(U_aphex{x, y});
		}
		return w;
	};

	EST_("real std::log")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += log(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2, M_car=1>")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1, 1>{}.template method<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=2>")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1>{}.template method<2>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=1>")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1>{}.template method<1>(x);
		}
		return w;
	};
	EST_("real logarithm... <N_lim=0>")
	{
		U_alpha w{};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f() + 1.0;
			w += logarithm_t< 1>{}.template method<0>(x);
		}
		return w;
	};
	EST_("real std::exp")
	{
		U_alpha o{1};
		for (int i{0x60}; ~--i;) {
			o *= exp(mt19937_f());
		}
		return o;
	
	};
	EST_("real antilogarithm... <N_lim=~0>")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			//\
			w *= logarithm_t<-1>{}.template method<~0>(mt19937_f());
			w *= exp(mt19937_f());
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2, M_car=1>")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1, 1>{}.template method<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=2>")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1>{}.template method<2>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=1>")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1>{}.template method<1>(x);
		}
		return w;
	
	};
	EST_("real antilogarithm... <N_lim=0>")
	{
		U_alpha w{1};
		for (int i{0x60}; ~--i;) {
			auto x = mt19937_f();
			w *= logarithm_t<-1>{}.template method<0>(x);
		}
		return w;
	
	};
};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
