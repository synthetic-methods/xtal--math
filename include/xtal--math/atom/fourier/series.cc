#pragma once
#include "./any.cc"

#include <unsupported/Eigen/FFT>



#include "./series.hh"
XTAL_ENV_(push)
namespace xtal::atom::math::fourier::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("fourier", "series")
{
	using U_fit = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	/**/
	TRY_("generation")
	{
		using W7 = series_t<U_alpha[0x7]>;
		using W8 = series_t<U_alpha[0x8]>;

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
	/**/
	TRY_("transformation")
	{
		U_sigma constexpr O = 1 << 3;
		U_sigma constexpr N = 1 << 3;
		U_sigma constexpr H = N >> 1;
		U_sigma constexpr M = N  - 1;

		using V_series = series_t<U_aphex[O]>;
		using U_series = series_t<U_aphex[N]>;
		V_series basis(constant_t<-1>{});

		U_series source;
		source[0] = source[M - 0] = U_aphex(0.0, 0.0);
		source[1] = source[M - 1] = U_aphex(1.0, 1.0);
		source[2] = source[M - 2] = U_aphex(3.0, 3.0);
		source[3] = source[M - 3] = U_aphex(4.0, 4.0);
		auto target = basis.transformation(source);

		std::vector<U_aphex> sourcery(source);
		std::vector<U_aphex> targetry(source);
		Eigen::FFT<U_alpha> fft; fft.fwd(targetry, sourcery);

		TRUE_(check_f<- 1>(target[0], targetry[0]));
		TRUE_(check_f<-12>(target[1], targetry[1]));
		TRUE_(check_f<- 1>(target[2], targetry[2]));
		TRUE_(check_f<-12>(target[3], targetry[3]));
		TRUE_(check_f<- 1>(target[4], targetry[4]));
		TRUE_(check_f<-12>(target[5], targetry[5]));
		TRUE_(check_f<- 1>(target[6], targetry[6]));
		TRUE_(check_f<-12>(target[7], targetry[7]));

		basis.template transform<-1>(target);
		TRUE_(check_f<- 1>(target[0], sourcery[0]));
		TRUE_(check_f<- 1>(target[1], sourcery[1]));
		TRUE_(check_f<- 1>(target[2], sourcery[2]));
		TRUE_(check_f<- 1>(target[3], sourcery[3]));
		TRUE_(check_f<- 1>(target[4], sourcery[4]));
		TRUE_(check_f<- 1>(target[5], sourcery[5]));
		TRUE_(check_f<- 1>(target[6], sourcery[6]));
		TRUE_(check_f<- 1>(target[7], sourcery[7]));
	}
	/***/
	/**/
	TRY_("convolution")
	{
		U_sigma constexpr N = 1 << 3;
		U_sigma constexpr M = N  - 1;

		using U_series = series_t<U_aphex[N]>;
		U_series basis(constant_t<-1>{});

		U_series lhs = {0, 1, 2, 0, 0, 0, 0, 0};
		U_series rhs = {1, 0, 1, 0, 0, 0, 0, 0};
		U_series xhs = {0, 1, 2, 1, 2, 0, 0, 0};
		U_series yhs = basis.convolution(lhs, rhs);
		_detail::apply_to<[] XTAL_1FN_(call) (bond::math::bit_trim_f<16>)>(yhs);
		TRUE_(xhs == yhs);

	}
	/***/
	/**/
	TRY_("multiplication")
	{
		using C4 = series_t<U_aphex[4]>;
		using D4 = series_t<U_aphex[4]>;
		
		TRUE_(C4{1000, 100, 10, 1} * C4{2000, 200, 20, 2} == C4{2000600, 400040, 60002, 8000});
		TRUE_(D4{1000, 100, 10, 1} * D4{2000, 200, 20, 2} == D4{2000600, 400040, 60002, 8000});

	}
	/***/
}
TAG_("fourier", "series", "fitness")
{
	using U_fit = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	auto mt19937_o = typename U_fit::MT19937{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (U_fit::mantissa_f(mt19937_o));

	U_sigma constexpr N = 1 << 8;
	U_sigma constexpr H = N >> 1;
	U_sigma constexpr M = N  - 1;

	using V_series = series_t<U_aphex[N]>;
	using U_series = series_t<U_aphex[N]>;
	V_series basis(constant_t<-1>{});

	U_series source;
	for (U_sigma i{}; i < H; ++i) {
		source[i] = source[M - i] = U_aphex(mt19937_f(), mt19937_f());
	}
	std::vector<U_aphex> sourcery(source);
	std::vector<U_aphex> targetry(source);
	Eigen::FFT<U_alpha> fft;
	/**/
	EST_("transform via Eigen/FFT")
	{
		fft.fwd(targetry, sourcery);// return targetry;

	};
	EST_("transform")
	{
		(void) basis.transform(source);// return &source;

	};
	EST_("transformation")
	{
		return basis.transformation(source);

	};
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
