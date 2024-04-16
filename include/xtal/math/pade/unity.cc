#pragma once
#include "./any.cc"
#include "./unity.hh"// testing...

#include "../dilate.hh"





XTAL_ENV_(push)
namespace xtal::math::pade::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_lim=0, int N_trim=0>
XTAL_FN2 unity__check_f(auto const &t)
XTAL_0EX
{
	int constexpr N_inf = -1;
	auto const u = unity_t<1>::template function<N_lim>(t);
	auto const v = unity_t<1>::template function<N_inf>(t);
	return bond::computrim_f<N_trim>(u) == bond::computrim_f<N_trim>(v);
}


////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("unity")
{
	using re = bond::realized;

	using T_sigma = typename re::sigma_t;
	using T_delta = typename re::delta_t;
	using T_alpha = typename re::alpha_t;
	using T_aphex = typename re::aphex_t;
	XTAL_LET_(T_alpha) two =  2;
	XTAL_LET_(T_alpha) ten = 10;

	using U_phi = algebra::differential::modular_t<T_alpha[2]>;

	auto mt19937_f = typename re::mt19937_t();
	mt19937_f.seed(Catch::rngSeed());

	TRY_("evaluation (floating-point)")
	{
		double const t0 = 1.125;

		TRUE_(unity__check_f<0,  2>(t0));
		TRUE_(unity__check_f<1,  6>(t0));
		TRUE_(unity__check_f<2, 14>(t0));
		TRUE_(unity__check_f<3, 20>(t0));
		TRUE_(unity__check_f<4, 28>(t0));
		TRUE_(unity__check_f<5, 37>(t0));
		TRUE_(unity__check_f<6, 47>(t0));
		TRUE_(unity__check_f<7, 49>(t0));

		EST_("evaluation <N_lim=-1>")
		{
			T_alpha t = 0.618, _t = 0.414;
			T_aphex z {1, 0};

			for (T_sigma i = 0x100; ~--i;) {
				z *= unity_t<1>::template function<-1>(t);
				t += _t;
				t -= _std::round(t);
			}
			return z;

		};
		EST_("bench <N_lim=4>")
		{
			T_alpha t = 0.618, _t = 0.414;
			T_aphex z {1, 0};

			for (T_sigma i = 0x100; ~--i;) {
				t += _t;
				t -= _std::round(t);
				z *= unity_t<1>::template function<4>(t);
			}
			return z;

		};
	}
	TRY_("evaluation (fixed-point)")
	{
		double const t0 = algebra::differential::modular_t<T_alpha[2]> {1.125, 0};

		TRUE_(unity__check_f<0,  2>(t0));
		TRUE_(unity__check_f<1,  6>(t0));
		TRUE_(unity__check_f<2, 14>(t0));
		TRUE_(unity__check_f<3, 20>(t0));
		TRUE_(unity__check_f<4, 28>(t0));
		TRUE_(unity__check_f<5, 37>(t0));
		TRUE_(unity__check_f<6, 47>(t0));
		TRUE_(unity__check_f<7, 49>(t0));

		EST_("bench <N_lim=4>")
		{
			U_phi t {0.618, 0.414};
			T_aphex z {1, 0};

			for (T_sigma i = 0x100; ~--i; ++t) {
				z *= unity_t<1>::template function<4>(t);
			}
			return z;

		};
	}
	TRY_("dilution")
	{
		double const t1 = 1.11;
		double const t2 = 2.22;

		int constexpr N_lim = 4;

		auto z = unity_t<1>::template function<N_lim>(t1);

		TRUE_(check_f<-1>(z, unity_t<1           >::template function<N_lim>(t1)));
		TRUE_(check_f<-1>(z, unity_t<1, dilute<1>>::template function<N_lim>(t2)));
		
		TRUE_(check_f<-1>(z, process::chain_t<unity<1>, dilute<1>>::template function<N_lim>(t2)));
		TRUE_(check_f<-1>(z, process::chain_t<unity<1>           >::template function<N_lim>(t1)));


	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
