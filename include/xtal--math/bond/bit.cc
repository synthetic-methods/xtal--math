#pragma once
#include "./any.cc"





#include "./bit.hh"
XTAL_ENV_(push)
namespace xtal::bond::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("bit", "constants")
{
	using T_fit = fit<unsigned>;
	using T_sigma = typename T_fit::sigma_type;
	using T_delta = typename T_fit::delta_type;

	TRY_("constant evaluation and access")
	{
		size_type constexpr N = 8;
		size_type constexpr M = bit_mask_f<N>();

		TRUE_(bit_prod_f<N>(M*1 - 2) == (M*1 - 2)%(M));
		TRUE_(bit_prod_f<N>(M*1 - 1) == (M*1 - 1)%(M));
		TRUE_(bit_prod_f<N>(M*1 + 0) == (M*1 + 0)%(M));
		TRUE_(bit_prod_f<N>(M*1 + 1) == (M*1 + 1)%(M));
		TRUE_(bit_prod_f<N>(M*1 + 2) == (M*1 + 2)%(M));

		TRUE_(bit_prod_f<N>(M*2 - 2) == (M*2 - 2)%(M));
		TRUE_(bit_prod_f<N>(M*2 - 1) == (M*2 - 1)%(M));
		TRUE_(bit_prod_f<N>(M*2 + 0) == (M*2 + 0)%(M));
		TRUE_(bit_prod_f<N>(M*2 + 1) == (M*2 + 1)%(M));
		TRUE_(bit_prod_f<N>(M*2 + 2) == (M*2 + 2)%(M));

		TRUE_(0xAAAAAAAAU ==  (~unsigned{}/3<<1));
		TRUE_(0x7FFFFFFFU ==  (~unsigned{}>>1));
		TRUE_(0x80000000U == ~(~unsigned{}>>1));

		TRUE_(0xAAAAAAAAU ==  T_fit::full.mask/3<<1);
		TRUE_(0x7FFFFFFFU == ~T_fit::sign.mask);
		TRUE_(0x80000000U ==  T_fit::sign.mask);

	}
}


////////////////////////////////////////////////////////////////////////////////

TAG_("bit", "bit_swap_f")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	TRY_("bit_swap_f<+1>( ordinal_q)")
	{
		T_delta x, y;

		x =  1; y =  2; TRUE_(bit_swap_f<+1>(x, y) ==  0); TRUE_(x ==  1); TRUE_(y ==  2);
		x =  2; y =  1; TRUE_(bit_swap_f<+1>(x, y) == -1); TRUE_(x ==  1); TRUE_(y ==  2);
		x = -1; y =  1; TRUE_(bit_swap_f<+1>(x, y) ==  0); TRUE_(x == -1); TRUE_(y ==  1);
		x =  1; y = -1; TRUE_(bit_swap_f<+1>(x, y) == -1); TRUE_(x == -1); TRUE_(y ==  1);

	}
	TRY_("bit_swap_f<+1>(cardinal_q)")
	{
		T_sigma x, y;

		x =  1; y =  2; TRUE_(bit_swap_f<+1>(x, y) ==  0); TRUE_(x ==  1); TRUE_(y ==  2);
		x =  2; y =  1; TRUE_(bit_swap_f<+1>(x, y) == -1); TRUE_(x ==  1); TRUE_(y ==  2);
	//	x = -1; y =  1; TRUE_(bit_swap_f<+1>(x, y) == -1); TRUE_(x == -1); TRUE_(y ==  1);
	//	x =  1; y = -1; TRUE_(bit_swap_f<+1>(x, y) ==  1); TRUE_(x ==  1); TRUE_(y == -1);

	}
}


////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("bit", "bit_zoom_f")
{
	using U_fit    = fit<>;
	using U_sigma  = typename U_fit::sigma_type;
	using U_delta  = typename U_fit::delta_type;
	using U_alpha  = typename U_fit::alpha_type;
	using U_aphex  = typename U_fit::aphex_type;
	auto mt19937_f = typename U_fit::MT19937(Catch::rngSeed());

	TRY_("bit_zoom_f evaluation")
	{
		U_alpha x{1'000};
		U_alpha y{0.001};
		U_aphex z{x, y};
		echo_(z*bit_zoom_f(z));

	}
}
/***/

////////////////////////////////////////////////////////////////////////////////

TAG_("bit", "bit_shift_f")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	TRY_("bit_shift_f evaluation")
	{
		TRUE_(444 == bit_shift_f(222,  1));
		TRUE_(111 == bit_shift_f(222, -1));

	}
}


////////////////////////////////////////////////////////////////////////////////

TAG_("bit", "bit_extremal_f")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	TRY_("bit_extremal_f(cardinal_q)")
	{
		T_sigma x, y;

		x = 15; y = 24; TRUE_(bit_extremal_f<-1>(x, y) ==   3);
		x = 42; y = 56; TRUE_(bit_extremal_f<-1>(x, y) ==  14);

		x = 15; y = 24; TRUE_(bit_extremal_f<+1>(x, y) == 120);
		x = 42; y = 56; TRUE_(bit_extremal_f<+1>(x, y) == 168);

	}
}


////////////////////////////////////////////////////////////////////////////////

TAG_("bit", "bit_trim_f")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	TRY_("trim")
	{
		TRUE_(bit_trim_f<2>(fit<>::patio_f(2, 1)) == 6.25);
		TRUE_(bit_trim_f<3>(fit<>::patio_f(1, 1)) == 3.125);
		TRUE_(bit_trim_f<4>(fit<>::patio_f(1, 2)) == 1.5625);

	}
}
/**/
TAG_("bit", "bit_reverse")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	TRY_("32:03")
	{
		using T_fit = fit<uint32_t>;
		T_sigma constexpr question = 0b011;
		T_sigma constexpr answer   = 0b110;
		T_sigma constexpr guess    = bit_reverse_f< 3>(question);
		TRUE_(answer == guess);

	}
	TRY_("16:16")
	{
		using T_fit = fit<uint16_t>;
		T_sigma constexpr question = 0b0100100011100101;
		T_sigma constexpr answer   = 0b1010011100010010;
		T_sigma constexpr guess    = bit_reverse_f<16>(question);
		TRUE_(answer == guess);

	}
	TRY_("16:12")
	{
		using T_fit = fit<uint16_t>;
		T_sigma constexpr question = 0b010010001110;
		T_sigma constexpr answer   = 0b011100010010;
		T_sigma constexpr guess    = bit_reverse_f<12>(question);
		TRUE_(answer == guess);

	}
	TRY_("8:8")
	{
		using T_fit = fit<uint8_t>;
		T_sigma constexpr question = 0b01001101;
		T_sigma constexpr answer   = 0b10110010;
		T_sigma constexpr guess    = bit_reverse_f< 8>(question);
		TRUE_(answer == guess);

	}
	TRY_("8:6")
	{
		using T_fit = fit<uint8_t>;
		T_sigma constexpr question = 0b010011;
		T_sigma constexpr answer   = 0b110010;
		T_sigma constexpr guess    = bit_reverse_f< 6>(question);
		TRUE_(answer == guess);

	}
}


////////////////////////////////////////////////////////////////////////////////

TAG_("bit_floor_f")
{
	using T_fit    = fit<>;
	using T_delta  = typename T_fit::delta_type;
	using T_sigma  = typename T_fit::sigma_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	TRY_("trial")
	{
		TRUE_(bit_depth_f(T_sigma(0b0000)) == 0/*0b0000*/);
		TRUE_(bit_depth_f(T_sigma(0b0001)) == 1/*0b0001*/);
		TRUE_(bit_depth_f(T_sigma(0b0010)) == 2/*0b0010*/);
		TRUE_(bit_depth_f(T_sigma(0b0011)) == 3/*0b0100*/);
		TRUE_(bit_depth_f(T_sigma(0b0100)) == 3/*0b0100*/);
		TRUE_(bit_depth_f(T_sigma(0b0101)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_sigma(0b0110)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_sigma(0b0111)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_sigma(0b1000)) == 4/*0b1000*/);

		TRUE_(bit_depth_f(T_delta(0b0000)) == 0/*0b0000*/);
		TRUE_(bit_depth_f(T_delta(0b0001)) == 1/*0b0001*/);
		TRUE_(bit_depth_f(T_delta(0b0010)) == 2/*0b0010*/);
		TRUE_(bit_depth_f(T_delta(0b0011)) == 3/*0b0100*/);
		TRUE_(bit_depth_f(T_delta(0b0100)) == 3/*0b0100*/);
		TRUE_(bit_depth_f(T_delta(0b0101)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_delta(0b0110)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_delta(0b0111)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_delta(0b1000)) == 4/*0b1000*/);

		TRUE_(bit_depth_f(T_alpha(0b0000)) == 0/*0b0000*/);
		TRUE_(bit_depth_f(T_alpha(0b0001)) == 1/*0b0001*/);
		TRUE_(bit_depth_f(T_alpha(0b0010)) == 2/*0b0010*/);
		TRUE_(bit_depth_f(T_alpha(0b0011)) == 3/*0b0100*/);
		TRUE_(bit_depth_f(T_alpha(0b0100)) == 3/*0b0100*/);
		TRUE_(bit_depth_f(T_alpha(0b0101)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_alpha(0b0110)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_alpha(0b0111)) == 4/*0b1000*/);
		TRUE_(bit_depth_f(T_alpha(0b1000)) == 4/*0b1000*/);

	//	TRUE_( 1 == bit_depth_f(typename T_fit::aphex_type{0, 1}));
	//	TRUE_( 1 == bit_depth_f(typename T_fit::aphex_type{0, 2}));
	//	TRUE_( 2 == bit_depth_f(typename T_fit::aphex_type{0, 3}));
	//	TRUE_( 2 == bit_depth_f(typename T_fit::aphex_type{0, 4}));
	//	TRUE_( 3 == bit_depth_f(typename T_fit::aphex_type{0, 5}));
	//	TRUE_( 3 == bit_depth_f(typename T_fit::aphex_type{0, 6}));
	//	TRUE_( 3 == bit_depth_f(typename T_fit::aphex_type{0, 7}));
	//	TRUE_( 3 == bit_depth_f(typename T_fit::aphex_type{0, 8}));

	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 0.0}) ==  0);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 0.5}) ==  0);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 1.0}) ==  0);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 1.5}) ==  1);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 2.0}) ==  1);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 2.5}) ==  2);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 3.0}) ==  2);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 3.5}) ==  2);
	//	TRUE_(bit_depth_f<0>(T_aphex{0.0, 4.0}) ==  2);

		TRUE_(bit_depth_f(T_aphex{3.0, 3.0}) ==  3);
	//	echo_(bit_depth_f(T_aphex{0.3, 0.3}), -1);//FIX:	Improper rounding?

		TRUE_(check_f<-1>(bit_exponent_f(T_alpha{0.5}), bit_exponent_f(T_aphex{0, 0.5})));
		TRUE_(check_f<-1>(bit_exponent_f(T_alpha{1.5}), bit_exponent_f(T_aphex{0, 1.5})));
		TRUE_(check_f<-1>(bit_exponent_f(T_alpha{2.5}), bit_exponent_f(T_aphex{0, 2.5})));
		TRUE_(check_f<-1>(bit_exponent_f(T_alpha{3.5}), bit_exponent_f(T_aphex{0, 3.5})));

	};
}


////////////////////////////////////////////////////////////////////////////////

TAG_("bit_separation_f")
{
	using T_fit    = fit<>;
	using T_delta  = typename T_fit::delta_type;
	using T_sigma  = typename T_fit::sigma_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	TRY_("trial")
	{
		TRUE_( 2.25 == bit_reparation_f(bit_separation_f( 2.25)));
		TRUE_( 1.25 == bit_reparation_f(bit_separation_f( 1.25)));
		TRUE_( 0.25 == bit_reparation_f(bit_separation_f( 0.25)));
		TRUE_( 0.75 == bit_reparation_f(bit_separation_f( 0.75)));
		TRUE_( 1.75 == bit_reparation_f(bit_separation_f( 1.75)));
		TRUE_( 2.75 == bit_reparation_f(bit_separation_f( 2.75)));

		TRUE_(-2.25 == bit_reparation_f(bit_separation_f(-2.25)));
		TRUE_(-1.25 == bit_reparation_f(bit_separation_f(-1.25)));
		TRUE_(-0.25 == bit_reparation_f(bit_separation_f(-0.25)));
		TRUE_(-0.75 == bit_reparation_f(bit_separation_f(-0.75)));
		TRUE_(-1.75 == bit_reparation_f(bit_separation_f(-1.75)));
		TRUE_(-2.75 == bit_reparation_f(bit_separation_f(-2.75)));

	};
}


////////////////////////////////////////////////////////////////////////////////

TAG_("semifractional")
{
	using T_fit   = fit<>;
	using T_delta = typename T_fit::delta_type;
	using T_sigma = typename T_fit::sigma_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;

	using U_fit   = fit<float>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	TRY_("trial")
	{
		using _qp = bond::fit<int>;
		
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-00), (+123456789.0e-00) - round(+123456789.0e-00)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-01), (+123456789.0e-01) - round(+123456789.0e-01)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-02), (+123456789.0e-02) - round(+123456789.0e-02)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-03), (+123456789.0e-03) - round(+123456789.0e-03)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-04), (+123456789.0e-04) - round(+123456789.0e-04)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-05), (+123456789.0e-05) - round(+123456789.0e-05)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-06), (+123456789.0e-06) - round(+123456789.0e-06)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-07), (+123456789.0e-07) - round(+123456789.0e-07)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-08), (+123456789.0e-08) - round(+123456789.0e-08)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-09), (+123456789.0e-09) - round(+123456789.0e-09)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-10), (+123456789.0e-10) - round(+123456789.0e-10)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-11), (+123456789.0e-11) - round(+123456789.0e-11)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-12), (+123456789.0e-12) - round(+123456789.0e-12)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-13), (+123456789.0e-13) - round(+123456789.0e-13)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-14), (+123456789.0e-14) - round(+123456789.0e-14)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-15), (+123456789.0e-15) - round(+123456789.0e-15)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-16), (+123456789.0e-16) - round(+123456789.0e-16)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-17), (+123456789.0e-17) - round(+123456789.0e-17)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-18), (+123456789.0e-18) - round(+123456789.0e-18)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+123456789.0e-19), (+123456789.0e-19) - round(+123456789.0e-19)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(+000000000.0e-09), (+000000000.0e-09) - round(+000000000.0e-09)));

		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-00), (-123456789.0e-00) - round(-123456789.0e-00)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-01), (-123456789.0e-01) - round(-123456789.0e-01)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-02), (-123456789.0e-02) - round(-123456789.0e-02)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-03), (-123456789.0e-03) - round(-123456789.0e-03)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-04), (-123456789.0e-04) - round(-123456789.0e-04)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-05), (-123456789.0e-05) - round(-123456789.0e-05)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-06), (-123456789.0e-06) - round(-123456789.0e-06)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-07), (-123456789.0e-07) - round(-123456789.0e-07)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-08), (-123456789.0e-08) - round(-123456789.0e-08)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-09), (-123456789.0e-09) - round(-123456789.0e-09)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-10), (-123456789.0e-10) - round(-123456789.0e-10)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-11), (-123456789.0e-11) - round(-123456789.0e-11)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-12), (-123456789.0e-12) - round(-123456789.0e-12)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-13), (-123456789.0e-13) - round(-123456789.0e-13)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-14), (-123456789.0e-14) - round(-123456789.0e-14)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-15), (-123456789.0e-15) - round(-123456789.0e-15)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-16), (-123456789.0e-16) - round(-123456789.0e-16)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-17), (-123456789.0e-17) - round(-123456789.0e-17)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-18), (-123456789.0e-18) - round(-123456789.0e-18)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-123456789.0e-19), (-123456789.0e-19) - round(-123456789.0e-19)));
		TRUE_(check_f<-26>(bit_fraction_f<U_alpha>()*bit_fraction_f<U_delta>(-000000000.0e-09), (-000000000.0e-09) - round(-000000000.0e-09)));

	};
}


////////////////////////////////////////////////////////////////////////////////

TAG_("fraction")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::MT19937(Catch::rngSeed());

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	TRY_("comparing implementations")
	{
		for (T_sigma i = 0x100; ~--i;) {
			T_alpha const u = ten*T_fit::mantissa_f(mt19937_f);
			TRUE_(check_f<16>(bit_fraction_f<T_alpha>(u), u - std::round(u)));
		}
	};
}
TAG_("bit trials")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_o = typename T_fit::MT19937{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (T_fit::mantissa_f(mt19937_o));

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	EST_("bond::bit_fraction_f\n   # - Round@#&\n   (*fixed-point*)")
	{
		T_delta w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto const u = ten*mt19937_f();
			w ^= bit_fraction_f<T_delta>(u);
		}
		return w;
	};
	EST_("bond::bit_fraction_f\n   # - Round@#&\n   (*floating-point*)")
	{
		T_delta w{};
		for (T_sigma i = 0x100; ~--i;) {
			auto const u = ten*mt19937_f();
			w ^= static_cast<T_delta>(T_fit::diplo_f(T_fit::full.depth)*(u - std::round(u)));
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
