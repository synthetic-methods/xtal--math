#pragma once
#include "./any.cc"





#include "./prime.hh"// testing...
XTAL_ENV_(push)
namespace xtal::bond::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("prime")
{
	using T_fit    = fit<>;
	using T_sigma  = typename T_fit::sigma_type;
	using T_delta  = typename T_fit::delta_type;
	using T_alpha  = typename T_fit::alpha_type;
	using T_aphex  = typename T_fit::aphex_type;
	auto mt19937_f = typename T_fit::mt19937_t(Catch::rngSeed());

	TRY_("prime_encode_f")
	{
		unsigned n_product = 2*3*3*3*7*13*13;
		TRUE_(prime_recode_f<0>(n_product) == 1);
		TRUE_(prime_recode_f<1>(n_product) == 3);
		TRUE_(prime_recode_f<2>(n_product) == 0);
		TRUE_(prime_recode_f<3>(n_product) == 1);
		TRUE_(prime_recode_f<4>(n_product) == 0);
		TRUE_(prime_recode_f<5>(n_product) == 2);
		TRUE_(prime_recode_f<6>(n_product) == 0);

	}
	TRY_("prime_recode_f")
	{
		unsigned n_product = 2*3*3*3*7*13*13;
		TRUE_(prime_recode_f<0>(n_product) == 1);
		TRUE_(prime_recode_f<1>(n_product) == 3);
		TRUE_(prime_recode_f<2>(n_product) == 0);
		TRUE_(prime_recode_f<3>(n_product) == 1);
		TRUE_(prime_recode_f<4>(n_product) == 0);
		TRUE_(prime_recode_f<5>(n_product) == 2);
		TRUE_(prime_recode_f<6>(n_product) == 0);

	}
	TRY_("prime_encode_f(i_primes) ")
	{
		TRUE_(prime_encode_f(0, 1, 0, 2) == 147);

		auto constexpr N_22335577   = prime_encode_f<2, 4>();
		auto constexpr N_2222333557 = prime_encode_f(bond::seek_reverse_s<4 + (1)>{});
		TRUE_(N_22335577   == 2*2*3*3*5*5*7*7    );
		TRUE_(N_2222333557 == 2*2*2*2*3*3*3*5*5*7);

		auto constexpr y0 = prime_encode_f<N_2222333557>(0, 2, 0,-1);
		TRUE_(check_f<8>(y0, 9./7.));

	//	TODO:
	//	-	Accommodate sign?
	//	-	Rework as `atom::primer`/`atom::math::primer`:
	//		-	Allow for `couple`/`dot` initialization.
	//		-	Provide operators for bag/set union/intersection/difference.

	//	using  primer = xtal::atom::math::primer_t<bond::bond::seek_reverse_s<4 + (1)>>;
	//	return primer::encode_f(0, 2, 0,-1);

	//	using  primer = xtal::atom::math::primer_t<bond::bond::seek_reverse_s<4 + (1)>>;
	//	primer primed(0, 2, 0,-1);
	//	primer primed(9./7.);

	//	i = primed.template decode<0>();

		//\
		auto x0 = y0;
		auto x0 = (extent_type) round(N_2222333557*y0);

		TRUE_(prime_recode_f<0, N_2222333557>(x0) ==  0);
		TRUE_(prime_recode_f<1, N_2222333557>(x0) ==  2);
		TRUE_(prime_recode_f<2, N_2222333557>(x0) ==  0);
		TRUE_(prime_recode_f<3, N_2222333557>(x0) == -1);
		TRUE_(prime_recode_f<4, N_2222333557>(x0) ==  0);
		TRUE_(prime_recode_f<5, N_2222333557>(x0) ==  0);
		TRUE_(prime_recode_f<6, N_2222333557>(x0) ==  0);
	}
	TRY_("prime_f<N_index>(n_power) ")
	{
		TRUE_(prime_f<  0>( 2) ==       4);
		TRUE_(prime_f<  1>( 2) ==       9);
		TRUE_(prime_f<  2>( 2) ==      25);

		TRUE_(prime_f<  0>(10) ==    1024);
		TRUE_(prime_f<  1>(10) ==   59049);
		TRUE_(prime_f<  2>(10) == 9765625);

	}
	TRY_("prime_f<N_index>()")
	{
		TRUE_(prime_f<  0>() ==   2);
		TRUE_(prime_f<  1>() ==   3);
		TRUE_(prime_f<  2>() ==   5);
		TRUE_(prime_f<  3>() ==   7);
		TRUE_(prime_f<  4>() ==  11);
		TRUE_(prime_f<  5>() ==  13);
		TRUE_(prime_f<  6>() ==  17);
		TRUE_(prime_f<  7>() ==  19);
		TRUE_(prime_f<  8>() ==  23);
		TRUE_(prime_f<  9>() ==  29);
		TRUE_(prime_f< 10>() ==  31);
		TRUE_(prime_f< 11>() ==  37);
		TRUE_(prime_f< 12>() ==  41);
		TRUE_(prime_f< 13>() ==  43);
		TRUE_(prime_f< 14>() ==  47);
		TRUE_(prime_f< 15>() ==  53);
		TRUE_(prime_f< 16>() ==  59);
		TRUE_(prime_f< 17>() ==  61);
		TRUE_(prime_f< 18>() ==  67);
		TRUE_(prime_f< 19>() ==  71);
		TRUE_(prime_f< 20>() ==  73);
		TRUE_(prime_f< 21>() ==  79);
		TRUE_(prime_f< 22>() ==  83);
		TRUE_(prime_f< 23>() ==  89);
		TRUE_(prime_f< 24>() ==  97);
		TRUE_(prime_f< 25>() == 101);
		TRUE_(prime_f< 26>() == 103);
		TRUE_(prime_f< 27>() == 107);
		TRUE_(prime_f< 28>() == 109);
		TRUE_(prime_f< 29>() == 113);
		TRUE_(prime_f< 30>() == 127);
		TRUE_(prime_f< 31>() == 131);
		TRUE_(prime_f< 32>() == 137);
		TRUE_(prime_f< 33>() == 139);
		TRUE_(prime_f< 34>() == 149);
		TRUE_(prime_f< 35>() == 151);
		TRUE_(prime_f< 36>() == 157);
		TRUE_(prime_f< 37>() == 163);
		TRUE_(prime_f< 38>() == 167);
		TRUE_(prime_f< 39>() == 173);
		TRUE_(prime_f< 40>() == 179);
		TRUE_(prime_f< 41>() == 181);
		TRUE_(prime_f< 42>() == 191);
		TRUE_(prime_f< 43>() == 193);
		TRUE_(prime_f< 44>() == 197);
		TRUE_(prime_f< 45>() == 199);
		TRUE_(prime_f< 46>() == 211);
		TRUE_(prime_f< 47>() == 223);
		TRUE_(prime_f< 48>() == 227);
		TRUE_(prime_f< 49>() == 229);
		TRUE_(prime_f< 50>() == 233);
		TRUE_(prime_f< 51>() == 239);
		TRUE_(prime_f< 52>() == 241);
		TRUE_(prime_f< 53>() == 251);
		TRUE_(prime_f< 54>() == 257);
		TRUE_(prime_f< 55>() == 263);
		TRUE_(prime_f< 56>() == 269);
		TRUE_(prime_f< 57>() == 271);
		TRUE_(prime_f< 58>() == 277);
		TRUE_(prime_f< 59>() == 281);
		TRUE_(prime_f< 60>() == 283);
		TRUE_(prime_f< 61>() == 293);
		TRUE_(prime_f< 62>() == 307);
		TRUE_(prime_f< 63>() == 311);
		TRUE_(prime_f< 64>() == 313);
		TRUE_(prime_f< 65>() == 317);
		TRUE_(prime_f< 66>() == 331);
		TRUE_(prime_f< 67>() == 337);
		TRUE_(prime_f< 68>() == 347);
		TRUE_(prime_f< 69>() == 349);
		TRUE_(prime_f< 70>() == 353);
		TRUE_(prime_f< 71>() == 359);
		TRUE_(prime_f< 72>() == 367);
		TRUE_(prime_f< 73>() == 373);
		TRUE_(prime_f< 74>() == 379);
		TRUE_(prime_f< 75>() == 383);
		TRUE_(prime_f< 76>() == 389);
		TRUE_(prime_f< 77>() == 397);
		TRUE_(prime_f< 78>() == 401);
		TRUE_(prime_f< 79>() == 409);
		TRUE_(prime_f< 80>() == 419);
		TRUE_(prime_f< 81>() == 421);
		TRUE_(prime_f< 82>() == 431);
		TRUE_(prime_f< 83>() == 433);
		TRUE_(prime_f< 84>() == 439);
		TRUE_(prime_f< 85>() == 443);
		TRUE_(prime_f< 86>() == 449);
		TRUE_(prime_f< 87>() == 457);
		TRUE_(prime_f< 88>() == 461);
		TRUE_(prime_f< 89>() == 463);
		TRUE_(prime_f< 90>() == 467);
		TRUE_(prime_f< 91>() == 479);
		TRUE_(prime_f< 92>() == 487);
		TRUE_(prime_f< 93>() == 491);
		TRUE_(prime_f< 94>() == 499);
		TRUE_(prime_f< 95>() == 503);
		TRUE_(prime_f< 96>() == 509);
		TRUE_(prime_f< 97>() == 521);
		TRUE_(prime_f< 98>() == 523);
		TRUE_(prime_f< 99>() == 541);
		TRUE_(prime_f<100>() == 547);
		TRUE_(prime_f<101>() == 557);
		TRUE_(prime_f<102>() == 563);
		TRUE_(prime_f<103>() == 569);
		TRUE_(prime_f<104>() == 571);
		TRUE_(prime_f<105>() == 577);
		TRUE_(prime_f<106>() == 587);
		TRUE_(prime_f<107>() == 593);
		TRUE_(prime_f<108>() == 599);
		TRUE_(prime_f<109>() == 601);
		TRUE_(prime_f<110>() == 607);
		TRUE_(prime_f<111>() == 613);
		TRUE_(prime_f<112>() == 617);
		TRUE_(prime_f<113>() == 619);
		TRUE_(prime_f<114>() == 631);
		TRUE_(prime_f<115>() == 641);
		TRUE_(prime_f<116>() == 643);
		TRUE_(prime_f<117>() == 647);
		TRUE_(prime_f<118>() == 653);
		TRUE_(prime_f<119>() == 659);
		TRUE_(prime_f<120>() == 661);
		TRUE_(prime_f<121>() == 673);
		TRUE_(prime_f<122>() == 677);
		TRUE_(prime_f<123>() == 683);
		TRUE_(prime_f<124>() == 691);
		TRUE_(prime_f<125>() == 701);
		TRUE_(prime_f<126>() == 709);
		TRUE_(prime_f<127>() == 719);

	}
	TRY_("prime_root_f")
	{
		TRUE_(prime_root_f(prime_f< 0>(), prime_root_f(prime_f< 0>())));
		TRUE_(prime_root_f(prime_f< 1>(), prime_root_f(prime_f< 1>())));
		TRUE_(prime_root_f(prime_f< 2>(), prime_root_f(prime_f< 2>())));
		TRUE_(prime_root_f(prime_f< 3>(), prime_root_f(prime_f< 3>())));
		TRUE_(prime_root_f(prime_f< 4>(), prime_root_f(prime_f< 4>())));
		TRUE_(prime_root_f(prime_f< 5>(), prime_root_f(prime_f< 5>())));
		TRUE_(prime_root_f(prime_f< 6>(), prime_root_f(prime_f< 6>())));
		TRUE_(prime_root_f(prime_f< 7>(), prime_root_f(prime_f< 7>())));
		TRUE_(prime_root_f(prime_f< 8>(), prime_root_f(prime_f< 8>())));
		TRUE_(prime_root_f(prime_f< 9>(), prime_root_f(prime_f< 9>())));
		TRUE_(prime_root_f(prime_f<10>(), prime_root_f(prime_f<10>())));
		TRUE_(prime_root_f(prime_f<11>(), prime_root_f(prime_f<11>())));
		TRUE_(prime_root_f(prime_f<12>(), prime_root_f(prime_f<12>())));
		TRUE_(prime_root_f(prime_f<13>(), prime_root_f(prime_f<13>())));
		TRUE_(prime_root_f(prime_f<14>(), prime_root_f(prime_f<14>())));
		TRUE_(prime_root_f(prime_f<15>(), prime_root_f(prime_f<15>())));
		TRUE_(prime_root_f(prime_f<16>(), prime_root_f(prime_f<16>())));
		TRUE_(prime_root_f(prime_f<17>(), prime_root_f(prime_f<17>())));
		TRUE_(prime_root_f(prime_f<18>(), prime_root_f(prime_f<18>())));
		TRUE_(prime_root_f(prime_f<19>(), prime_root_f(prime_f<19>())));
		TRUE_(prime_root_f(prime_f<20>(), prime_root_f(prime_f<20>())));
		TRUE_(prime_root_f(prime_f<21>(), prime_root_f(prime_f<21>())));
		TRUE_(prime_root_f(prime_f<22>(), prime_root_f(prime_f<22>())));
		TRUE_(prime_root_f(prime_f<23>(), prime_root_f(prime_f<23>())));
		TRUE_(prime_root_f(prime_f<24>(), prime_root_f(prime_f<24>())));
		TRUE_(prime_root_f(prime_f<25>(), prime_root_f(prime_f<25>())));
		TRUE_(prime_root_f(prime_f<26>(), prime_root_f(prime_f<26>())));
		TRUE_(prime_root_f(prime_f<27>(), prime_root_f(prime_f<27>())));
		TRUE_(prime_root_f(prime_f<28>(), prime_root_f(prime_f<28>())));
		TRUE_(prime_root_f(prime_f<29>(), prime_root_f(prime_f<29>())));
		TRUE_(prime_root_f(prime_f<30>(), prime_root_f(prime_f<30>())));
		TRUE_(prime_root_f(prime_f<31>(), prime_root_f(prime_f<31>())));
		TRUE_(prime_root_f(prime_f<32>(), prime_root_f(prime_f<32>())));
		TRUE_(prime_root_f(prime_f<33>(), prime_root_f(prime_f<33>())));
		TRUE_(prime_root_f(prime_f<34>(), prime_root_f(prime_f<34>())));
		TRUE_(prime_root_f(prime_f<35>(), prime_root_f(prime_f<35>())));
		TRUE_(prime_root_f(prime_f<36>(), prime_root_f(prime_f<36>())));
		TRUE_(prime_root_f(prime_f<37>(), prime_root_f(prime_f<37>())));
		TRUE_(prime_root_f(prime_f<38>(), prime_root_f(prime_f<38>())));
		TRUE_(prime_root_f(prime_f<39>(), prime_root_f(prime_f<39>())));
		TRUE_(prime_root_f(prime_f<40>(), prime_root_f(prime_f<40>())));
		TRUE_(prime_root_f(prime_f<41>(), prime_root_f(prime_f<41>())));
		TRUE_(prime_root_f(prime_f<42>(), prime_root_f(prime_f<42>())));
		TRUE_(prime_root_f(prime_f<43>(), prime_root_f(prime_f<43>())));
		TRUE_(prime_root_f(prime_f<44>(), prime_root_f(prime_f<44>())));
		TRUE_(prime_root_f(prime_f<45>(), prime_root_f(prime_f<45>())));
		TRUE_(prime_root_f(prime_f<46>(), prime_root_f(prime_f<46>())));
		TRUE_(prime_root_f(prime_f<47>(), prime_root_f(prime_f<47>())));
		TRUE_(prime_root_f(prime_f<48>(), prime_root_f(prime_f<48>())));
		TRUE_(prime_root_f(prime_f<49>(), prime_root_f(prime_f<49>())));

	//	TRUE_(prime_root_f(prime_f<50>(), prime_root_f(prime_f<50>())));// `233`
	//	TRUE_(prime_root_f(prime_f<51>(), prime_root_f(prime_f<51>())));// `239`
	//	TRUE_(prime_root_f(prime_f<52>(), prime_root_f(prime_f<52>())));// `241`
	//	TRUE_(prime_root_f(prime_f<53>(), prime_root_f(prime_f<53>())));// `251`

		TRUE_(prime_root_f(prime_f<54>(), prime_root_f(prime_f<54>())));// `257`

	}
}


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
