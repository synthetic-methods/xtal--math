#pragma once
#include "./any.cc"
#include "./root.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("root")
{
	using _fit = bond::fit<>;
	using T_sigma  = typename _fit::sigma_type;
	using T_delta  = typename _fit::delta_type;
	using T_alpha  = typename _fit::alpha_type;
	using T_aphex  = typename _fit::aphex_type;

	auto constexpr N_half_root2 = root_f<-2>(_fit::diplo_1);
	auto constexpr N_half_root3 = root_f<-3>(_fit::diplo_1);
	auto constexpr N_half_root4 = root_f<-4>(_fit::diplo_1);
	auto constexpr N_half_root5 = root_f<-5>(_fit::diplo_1);

	TRY_("evaluation (constexpr)")
	{
	//	NOTE: To stop GCC from choking on `REQUIRE`... (possible issue with ` __attribute__((optimize(0))) constexpr`?)
		static_assert(check_f<- 1>(N_half_root2, T_alpha{0.7071067811865475244008443621048490L}));
		static_assert(check_f<- 1>(N_half_root3, T_alpha{0.7937005259840997373758528196361541L}));
		static_assert(check_f<- 1>(N_half_root4, T_alpha{0.8408964152537145430311254762332149L}));
		static_assert(check_f<- 1>(N_half_root5, T_alpha{0.8705505632961241391362700174797461L}));

	}
	TRY_("evaluation")
	{
		TRUE_(check_f<- 2>(root_f<-1, 1>(T_aphex{1<<24, 1<<24}), one/T_aphex{1<<24, 1<<24}));
		TRUE_(check_f<- 1>(root_f<-2>(T_aphex{1, 1}), one/root_f< 2>(T_aphex{1, 1})));

		TRUE_(check_f<-30>(pow(0.5, _fit::ratio_f(1,  3)), root_f< 3>(0.5)));
		TRUE_(check_f<-30>(pow(0.5, _fit::ratio_f(1, -3)), root_f<-3>(0.5)));

		TRUE_(check_f<-30>(pow(0.5, _fit::ratio_f(1,  5)), root_f< 5>(0.5)));
		TRUE_(check_f<-30>(pow(0.5, _fit::ratio_f(1, -5)), root_f<-5>(0.5)));

		TRUE_(check_f<-30>(pow(0.5, _fit::ratio_f(1,  7)), root_f< 7>(0.5)));
		TRUE_(check_f<-30>(pow(0.5, _fit::ratio_f(1, -7)), root_f<-7>(0.5)));

		TRUE_(check_f<-30>(1.0, cbrt(-0.500)*root_f<-3>(-0.500)));
		TRUE_(check_f<-30>(1.0, cbrt(-0.250)*root_f<-3>(-0.250)));
		TRUE_(check_f<-30>(1.0, cbrt(-0.125)*root_f<-3>(-0.125)));

		TRUE_(check_f<- 1>(root_f< 4>(0.5),     sqrt(sqrt(0.5))));
		TRUE_(check_f<- 1>(root_f<-4>(0.5), 1.0/sqrt(sqrt(0.5))));

		TRUE_(check_f<- 1>(root_f< 8>(0.5),     sqrt(sqrt(sqrt(0.5)))));
		TRUE_(check_f<- 1>(root_f<-8>(0.5), 1.0/sqrt(sqrt(sqrt(0.5)))));

		TRUE_(check_f<- 1>(root_f<-2>(pow(T_aphex {2, 3}, -2.0)), T_aphex {2, 3}));
		TRUE_(check_f<- 1>(pow(T_aphex {2, 3}, 0.5), root_f< 2>(T_aphex {2, 3})));

		TRUE_(check_f<- 1>(1.0/sqrt(2.2345268795805384), root_f<-2>(_std::complex{2.2345268795805384, 0.0}).real()));
	}
	TRY_("punctured evaluation")
	{
		TRUE_(0.0 == root_f<-1, 1>(0.0) - root_f<-1, 1>(0.0));
		TRUE_(0.0 == root_f<-2, 1>(0.0) - root_f<-2, 1>(0.0));

	}

}
TAG_("root trials")
{
	using _fit = bond::fit<>;
	using T_sigma  = typename _fit::sigma_type;
	using T_delta  = typename _fit::delta_type;
	using T_alpha  = typename _fit::alpha_type;
	using T_aphex  = typename _fit::aphex_type;
	auto mt19937_o = typename _fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (_fit::mantissa_f(mt19937_o));

	auto constexpr N_half_root2 = root_f<-2>(_fit::diplo_1);
	auto constexpr N_half_root3 = root_f<-3>(_fit::diplo_1);
	auto constexpr N_half_root4 = root_f<-4>(_fit::diplo_1);
	auto constexpr N_half_root5 = root_f<-5>(_fit::diplo_1);

	EST_("root<-3;  3>\n   #^(1/3)&\n   (*native real*)")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= one/cbrt(mt19937_f() + one);
		}
		return w;

	};
	EST_("root<-3;  3>\n   #^(1/3)&\n   (*approx real*)")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= root_t<-3>::template method_f<3>(mt19937_f() + one);
		}
		return w;

	};
	EST_("root<-3;  2>\n   #^(1/3)&\n   (*approx real*)")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= root_t<-3>::template method_f<2>(mt19937_f() + one);
		}
		return w;

	};
	EST_("root<-3;  1>\n   #^(1/3)&\n   (*approx real*)")
	{
		double w{1};
		for (int i = 0x100; ~--i;) {
			w *= root_t<-3>::template method_f<1>(mt19937_f() + one);
		}
		return w;

	};
};
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
