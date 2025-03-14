#pragma once
#include "./any.cc"
#include "./series.hh"// testing...





XTAL_ENV_(push)
namespace xtal::atom::math::fourier::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//atic_assert(_xtd::trivially_initializable<series_t<float[2]>>);// Why?
//static_assert(_xtd::trivially_destructible <series_t<float[2]>>);
//static_assert(_xtd::trivially_copyable     <series_t<float[2]>>);
//static_assert(_xtd::trivially_movable      <series_t<float[2]>>);
//atic_assert(                     atomic_q<series_t<float[2]>>);// Why?


////////////////////////////////////////////////////////////////////////////////

TAG_("solid", "series")
{
	using _fit = bond::fit<>;
	using T_delta = typename _fit::delta_type;
	using T_sigma = typename _fit::sigma_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	/**/
	TRY_("initialization")
	{
		T_sigma constexpr N = 1 << 3;
		using V_series = series_t<T_alpha[N]>;
		using U_series = series_t<T_alpha[N]>;

		U_series baz(2.0);
		V_series bar = reinterpret_cast<V_series &>(baz);
		V_series foo = {1<<0, 1<<1, 1<<2, 1<<3, 1<<4, 1<<5, 1<<6, 1<<7};
		TRUE_(equal_f(foo, bar));
		
	}
	/***/
	/*/
	TRY_("generation of single")
	{
		using W7 = series_t<T_alpha[0x7]>;
		using W8 = series_t<T_alpha[0x8]>;

		W7 w7; w7.template generate<-1>(2.0);
		W8 w8; w8.template generate<-1>(2.0);

		TRUE_((1 << 1) == get<0>(w7));
		TRUE_((1 << 2) == get<1>(w7));
		TRUE_((1 << 3) == get<2>(w7));
		TRUE_((1 << 4) == get<3>(w7));
		TRUE_((1 << 5) == get<4>(w7));
		TRUE_((1 << 6) == get<5>(w7));
		TRUE_((1 << 7) == get<6>(w7));
	//	TRUE_((1 << 8) == get<7>(w7));

		TRUE_((1 << 1) == get<0>(w8));
		TRUE_((1 << 2) == get<1>(w8));
		TRUE_((1 << 3) == get<2>(w8));
		TRUE_((1 << 4) == get<3>(w8));
		TRUE_((1 << 5) == get<4>(w8));
		TRUE_((1 << 6) == get<5>(w8));
		TRUE_((1 << 7) == get<6>(w8));
		TRUE_((1 << 8) == get<7>(w8));

	}
	/***/
	/*/
	TRY_("generation of couple")
	{
		using U_alpha = couple_t<T_alpha[1<<1]>;
		using U_aphex = couple_t<T_aphex[1<<1]>;
		using W_aphex = series_t<U_aphex[1<<4]>;

		W_aphex w_aphex; w_aphex.generate(T_aphex{0, 1}, T_alpha{2.0});
		TRUE_(get<0>(w_aphex) == U_aphex {{  1.00000,  0.00000}, {  1.00000,- 0.00000}});
		TRUE_(get<1>(w_aphex) == U_aphex {{  0.00000,  2.00000}, {  0.00000,- 0.50000}});
		TRUE_(get<2>(w_aphex) == U_aphex {{- 4.00000,  0.00000}, {- 0.25000,- 0.00000}});
		TRUE_(get<3>(w_aphex) == U_aphex {{- 0.00000,- 8.00000}, {- 0.00000,  0.12500}});
		TRUE_(get<4>(w_aphex) == U_aphex {{ 16.00000,- 0.00000}, {  0.06250,  0.00000}});
		TRUE_(get<5>(w_aphex) == U_aphex {{  0.00000, 32.00000}, {  0.00000,- 0.03125}});

	}
	/***/
	/*/
	TRY_("transformation")
	{
		T_sigma constexpr O = 1 << 5;
		T_sigma constexpr N = 1 << 3;
		T_sigma constexpr M = N  - 1;

		using V_series = series_t<T_aphex[O]>;
		using U_series = series_t<T_aphex[N]>;
		V_series basis(constant_t<-1>{});

		U_series source;
		source[0] = source[M - 0] = T_aphex(0.0, 0.0);
		source[1] = source[M - 1] = T_aphex(1.0, 1.0);
		source[2] = source[M - 2] = T_aphex(3.0, 3.0);
		source[3] = source[M - 3] = T_aphex(4.0, 4.0);

		auto target = basis.transformation(source);
		TRUE_(check_f<-6>(target[0], T_aphex{ 0.1600000000000000e+2,  0.1600000000000000e+2}));
		TRUE_(check_f<-6>(target[1], T_aphex{-0.4828427124746192e+1, -0.1165685424949238e+2}));
		TRUE_(check_f<-6>(target[2], T_aphex{ 0.0000000000000000e+0,  0.0000000000000000e+0}));
		TRUE_(check_f<-6>(target[3], T_aphex{-0.3431457505076203e+0,  0.8284271247461885e+0}));
		TRUE_(check_f<-6>(target[4], T_aphex{ 0.0000000000000000e+0,  0.0000000000000000e+0}));
		TRUE_(check_f<-6>(target[5], T_aphex{ 0.8284271247461912e+0, -0.3431457505076203e+0}));
		TRUE_(check_f<-6>(target[6], T_aphex{ 0.0000000000000000e+0,  0.0000000000000000e+0}));
		TRUE_(check_f<-6>(target[7], T_aphex{-0.1165685424949238e+2, -0.4828427124746188e+1}));

	}
	/***/
	/*/
	TRY_("convolution")
	{
		T_sigma constexpr N = 1 << 3;
		T_sigma constexpr M = N  - 1;

		using U_series = series_t<T_aphex[N]>;
		U_series basis(constant_t<-1>{});

		U_series lhs = {0, 1, 2, 0, 0, 0, 0, 0};
		U_series rhs = {1, 0, 1, 0, 0, 0, 0, 0};
		U_series xhs = {0, 1, 2, 1, 2, 0, 0, 0};
		U_series yhs = basis.convolution(lhs, rhs);
		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(yhs);
		TRUE_(xhs == yhs);

	}
	/***/
	/*/
	TRY_("multiplication")
	{
		using C4 = series_t<T_aphex[4]>;
		using D4 = series_t<T_aphex[4]>;
		
		TRUE_(C4{1000, 100, 10, 1} * C4{2000, 200, 20, 2} == C4{2000600, 400040, 60002, 8000});
		TRUE_(D4{1000, 100, 10, 1} * D4{2000, 200, 20, 2} == D4{2000600, 400040, 60002, 8000});

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
