#pragma once
#include "./any.cc"
#include "./bit.hh"// testing...





XTAL_ENV_(push)
namespace xtal::bond::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("bit", "bit_trim_f")
{
	using T_fit = fit<>;
	using T_sigma = typename T_fit::sigma_type;
	using T_delta = typename T_fit::delta_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::mt19937_t(Catch::rngSeed());

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
	using T_fit = fit<>;
	using T_sigma = typename T_fit::sigma_type;
	using T_delta = typename T_fit::delta_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::mt19937_t(Catch::rngSeed());


	TRY_("32:03")
	{
		using T_fit = fit<uint32_t>;
		T_sigma question = 0b011;
		T_sigma answer   = 0b110;

		TRUE_(answer == bit_reverse_f<3>(question));

	}
	TRY_("16:16")
	{
		using T_fit = fit<uint16_t>;
		T_sigma question = 0b0100100011100101;
		T_sigma answer   = 0b1010011100010010;

		TRUE_(answer == bit_reverse_f<16>(question));

	}
	TRY_("16:12")
	{
		using T_fit = fit<uint16_t>;
		T_sigma question = 0b010010001110;
		T_sigma answer   = 0b011100010010;

		TRUE_(answer == bit_reverse_f<12>(question));

	}
	TRY_("8:8")
	{
		using T_fit = fit<uint8_t>;
		T_sigma question = 0b01001101;
		T_sigma answer   = 0b10110010;

		TRUE_(answer == bit_reverse_f<8>(question));

	}
	TRY_("8:6")
	{
		using T_fit = fit<uint8_t>;
		T_sigma question = 0b010011;
		T_sigma answer   = 0b110010;

		TRUE_(answer == bit_reverse_f<6>(question));

	}
}


////////////////////////////////////////////////////////////////////////////////

TAG_("bit_floor_f")
{
	using T_fit = fit<>;
	using T_delta = typename T_fit::delta_type;
	using T_sigma = typename T_fit::sigma_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::mt19937_t(Catch::rngSeed());

	TRY_("trial")
	{
		TRUE_( 9 == bit_floor_f<-4>(T_fit::diplo_f(9)));
		TRUE_( 8 == bit_floor_f<-4>(T_fit::diplo_f(8)));
		TRUE_( 7 == bit_floor_f<-4>(T_fit::diplo_f(7)));
		TRUE_( 6 == bit_floor_f<-4>(T_fit::diplo_f(6)));
		TRUE_( 5 == bit_floor_f<-4>(T_fit::diplo_f(5)));
		TRUE_( 4 == bit_floor_f<-4>(T_fit::diplo_f(4)));
		TRUE_( 3 == bit_floor_f<-4>(T_fit::diplo_f(3)));
		TRUE_( 2 == bit_floor_f<-4>(T_fit::diplo_f(2)));
		TRUE_( 1 == bit_floor_f<-4>(T_fit::diplo_f(1)));
		TRUE_( 0 == bit_floor_f<-4>(T_fit::diplo_f(0)));
		TRUE_( 0 == bit_floor_f<-4>(T_fit::haplo_f(0)));
		TRUE_(-1 == bit_floor_f<-4>(T_fit::haplo_f(1)));
		TRUE_(-2 == bit_floor_f<-4>(T_fit::haplo_f(2)));
		TRUE_(-3 == bit_floor_f<-4>(T_fit::haplo_f(3)));
		TRUE_(-4 == bit_floor_f<-4>(T_fit::haplo_f(4)));
		TRUE_(-4 == bit_floor_f<-4>(T_fit::haplo_f(5)));
		TRUE_(-4 == bit_floor_f<-4>(T_fit::haplo_f(6)));
		TRUE_(-4 == bit_floor_f<-4>(T_fit::haplo_f(7)));
		TRUE_(-4 == bit_floor_f<-4>(T_fit::haplo_f(8)));
		TRUE_(-4 == bit_floor_f<-4>(T_fit::haplo_f(9)));
		TRUE_(-4 == bit_floor_f<-4>(T_fit::alpha_f(0)));


		TRUE_(0 == bit_floor_f<0>(0));
		TRUE_(0 == bit_floor_f<0>(1));
		TRUE_(1 == bit_floor_f<0>(2));
		TRUE_(0 == bit_floor_f<0>(0.0));
		TRUE_(0 == bit_floor_f<0>(1.0));
		TRUE_(1 == bit_floor_f<0>(2.0));
	//	TRUE_(bit_floor_f(0.0) == bit_floor_f(0));
		TRUE_(bit_floor_f(1.0) == bit_floor_f(1));
		TRUE_(bit_floor_f(2.0) == bit_floor_f(2));
		TRUE_(bit_floor_f(3.0) == bit_floor_f(3));
		TRUE_(bit_floor_f(4.0) == bit_floor_f(4));
		TRUE_(bit_floor_f(5.0) == bit_floor_f(5));
		TRUE_(bit_floor_f(6.0) == bit_floor_f(6));
		TRUE_(bit_floor_f(7.0) == bit_floor_f(7));
		TRUE_(bit_floor_f(8.0) == bit_floor_f(8));
		TRUE_(bit_floor_f(9.0) == bit_floor_f(9));

	//	TRUE_(bit_floor_f(0.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(0));
		TRUE_(bit_floor_f(1.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(1));
		TRUE_(bit_floor_f(2.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(2));
		TRUE_(bit_floor_f(3.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(3));
		TRUE_(bit_floor_f(4.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(4));
		TRUE_(bit_floor_f(5.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(5));
		TRUE_(bit_floor_f(6.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(6));
		TRUE_(bit_floor_f(7.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(7));
		TRUE_(bit_floor_f(8.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(8));
		TRUE_(bit_floor_f(9.0*T_fit::diplo_1*T_fit::dnsilon_f(1)) == bit_ceiling_f(9));

	//	TRUE_(bit_ceiling_f(0) == _std::bit_width(0U - 1));
		TRUE_(bit_ceiling_f(1) == _std::bit_width(1U - 1));
		TRUE_(bit_ceiling_f(2) == _std::bit_width(2U - 1));
		TRUE_(bit_ceiling_f(3) == _std::bit_width(3U - 1));
		TRUE_(bit_ceiling_f(4) == _std::bit_width(4U - 1));
		TRUE_(bit_ceiling_f(5) == _std::bit_width(5U - 1));
		TRUE_(bit_ceiling_f(6) == _std::bit_width(6U - 1));
		TRUE_(bit_ceiling_f(7) == _std::bit_width(7U - 1));
		TRUE_(bit_ceiling_f(8) == _std::bit_width(8U - 1));
		TRUE_(bit_ceiling_f(9) == _std::bit_width(9U - 1));

		TRUE_(-1 == bond::math::bit_floor_f<-2>(-0.50));
		TRUE_(-1 == bond::math::bit_floor_f<-2>(+0.50));
		TRUE_(-2 == bond::math::bit_floor_f<-2>(-0.25));
		TRUE_(-2 == bond::math::bit_floor_f<-2>(+0.25));
		TRUE_( 1 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 1}));
		TRUE_( 1 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 2}));
		TRUE_( 2 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 3}));
		TRUE_( 2 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 4}));
		TRUE_( 3 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 5}));
		TRUE_( 3 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 6}));
		TRUE_( 3 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 7}));
		TRUE_( 3 == bond::math::bit_ceiling_f(typename T_fit::aphex_type{0, 8}));

		TRUE_(bit_floor_f(0.1) == -4);
		TRUE_(bit_floor_f(0.2) == -3);
		TRUE_(bit_floor_f(0.3) == -2);
		TRUE_(bit_floor_f(0.4) == -2);
		TRUE_(bit_floor_f(0.5) == -1);
		TRUE_(bit_floor_f(0.6) == -1);
		TRUE_(bit_floor_f(0.7) == -1);
		TRUE_(bit_floor_f(0.8) == -1);
		TRUE_(bit_floor_f(0.9) == -1);
		TRUE_(bit_floor_f(1.0) ==  0);
		TRUE_(bit_floor_f(2.0) ==  1);
		TRUE_(bit_floor_f(3.0) ==  1);
		TRUE_(bit_floor_f(4.0) ==  2);
		TRUE_(bit_floor_f(5.0) ==  2);
		TRUE_(bit_floor_f(6.0) ==  2);
		TRUE_(bit_floor_f(7.0) ==  2);
		TRUE_(bit_floor_f(8.0) ==  3);
		TRUE_(bit_floor_f(9.0) ==  3);

		TRUE_(bit_floor_f  (T_aphex{3.0, 3.0}) ==  2);
		TRUE_(bit_floor_f  (T_aphex{0.3, 0.3}) == -2);
	
		TRUE_(bit_ceiling_f(T_aphex{3.0, 3.0}) ==  3);
	//	echo_(bit_ceiling_f(T_aphex{0.3, 0.3}), -1);//FIX:	Improper rounding?


	};
}


////////////////////////////////////////////////////////////////////////////////

TAG_("bit_representation_f")
{
	using T_fit = fit<>;
	using T_delta = typename T_fit::delta_type;
	using T_sigma = typename T_fit::sigma_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::mt19937_t(Catch::rngSeed());

	TRY_("trial")
	{
		TRUE_( 2.25 == bit_presentation_f(bit_representation_f( 2.25)));
		TRUE_( 1.25 == bit_presentation_f(bit_representation_f( 1.25)));
		TRUE_( 0.25 == bit_presentation_f(bit_representation_f( 0.25)));
		TRUE_( 0.75 == bit_presentation_f(bit_representation_f( 0.75)));
		TRUE_( 1.75 == bit_presentation_f(bit_representation_f( 1.75)));
		TRUE_( 2.75 == bit_presentation_f(bit_representation_f( 2.75)));

		TRUE_(-2.25 == bit_presentation_f(bit_representation_f(-2.25)));
		TRUE_(-1.25 == bit_presentation_f(bit_representation_f(-1.25)));
		TRUE_(-0.25 == bit_presentation_f(bit_representation_f(-0.25)));
		TRUE_(-0.75 == bit_presentation_f(bit_representation_f(-0.75)));
		TRUE_(-1.75 == bit_presentation_f(bit_representation_f(-1.75)));
		TRUE_(-2.75 == bit_presentation_f(bit_representation_f(-2.75)));

	};
}


////////////////////////////////////////////////////////////////////////////////

TAG_("semifractional")
{
	using T_fit = fit<>;
	using T_delta = typename T_fit::delta_type;
	using T_sigma = typename T_fit::sigma_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;

	using U_fit = fit<float>;
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
	using T_fit = fit<>;
	using T_sigma = typename T_fit::sigma_type;
	using T_delta = typename T_fit::delta_type;
	using T_alpha = typename T_fit::alpha_type;
	using T_aphex = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::mt19937_t(Catch::rngSeed());

	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	TRY_("comparing implementations")
	{
		for (T_sigma i = 0x100; ~--i;) {
			T_alpha const u = ten*T_fit::mantissa_f(mt19937_f);
			TRUE_(check_f<16>(bit_fraction_f<T_alpha>(u), u - _std::round(u)));
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
	auto mt19937_o = typename T_fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
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
			w ^= static_cast<T_delta>(T_fit::diplo_f(T_fit::full.depth)*(u - _std::round(u)));
		}
		return w;
	};
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
