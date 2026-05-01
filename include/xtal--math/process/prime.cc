#pragma once
#include "./any.cc"





#include "./prime.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("prime")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;


	TRY_("prime_f<_std::array<...>>")
	{
		using U_nom = _std::array<U_delta, 4>;
		U_nom constexpr M_nom{3, 3, 3, 3};
		U_nom constexpr   vec{1, 2, 0,-3};
		
		using code = prime_t<M_nom>;
		auto  constexpr val = code::method(vec);

		TRUE_(check_f<-1>(U_fit::ratio_f(2*3*3, 7*7*7), val));
		TRUE_(vec == code::method(val));

	}
	TRY_("prime_f<-1>(prime_f<1>(...))")
	{
		#pragma unroll
		for (int i{0}; i < 0x80; ++i) {
			TRUE_(i == prime_f<-1>(prime_f<1>(i)));
		}
	}
	TRY_("prime_f<1>(...)")
	{
		TRUE_(prime_f<1>(  0) ==   2);
		TRUE_(prime_f<1>(  1) ==   3);
		TRUE_(prime_f<1>(  2) ==   5);
		TRUE_(prime_f<1>(  3) ==   7);
		TRUE_(prime_f<1>(  4) ==  11);
		TRUE_(prime_f<1>(  5) ==  13);
		TRUE_(prime_f<1>(  6) ==  17);
		TRUE_(prime_f<1>(  7) ==  19);
		TRUE_(prime_f<1>(  8) ==  23);
		TRUE_(prime_f<1>(  9) ==  29);
		TRUE_(prime_f<1>( 10) ==  31);
		TRUE_(prime_f<1>( 11) ==  37);
		TRUE_(prime_f<1>( 12) ==  41);
		TRUE_(prime_f<1>( 13) ==  43);
		TRUE_(prime_f<1>( 14) ==  47);
		TRUE_(prime_f<1>( 15) ==  53);
		TRUE_(prime_f<1>( 16) ==  59);
		TRUE_(prime_f<1>( 17) ==  61);
		TRUE_(prime_f<1>( 18) ==  67);
		TRUE_(prime_f<1>( 19) ==  71);
		TRUE_(prime_f<1>( 20) ==  73);
		TRUE_(prime_f<1>( 21) ==  79);
		TRUE_(prime_f<1>( 22) ==  83);
		TRUE_(prime_f<1>( 23) ==  89);
		TRUE_(prime_f<1>( 24) ==  97);
		TRUE_(prime_f<1>( 25) == 101);
		TRUE_(prime_f<1>( 26) == 103);
		TRUE_(prime_f<1>( 27) == 107);
		TRUE_(prime_f<1>( 28) == 109);
		TRUE_(prime_f<1>( 29) == 113);
		TRUE_(prime_f<1>( 30) == 127);
		TRUE_(prime_f<1>( 31) == 131);
		TRUE_(prime_f<1>( 32) == 137);
		TRUE_(prime_f<1>( 33) == 139);
		TRUE_(prime_f<1>( 34) == 149);
		TRUE_(prime_f<1>( 35) == 151);
		TRUE_(prime_f<1>( 36) == 157);
		TRUE_(prime_f<1>( 37) == 163);
		TRUE_(prime_f<1>( 38) == 167);
		TRUE_(prime_f<1>( 39) == 173);
		TRUE_(prime_f<1>( 40) == 179);
		TRUE_(prime_f<1>( 41) == 181);
		TRUE_(prime_f<1>( 42) == 191);
		TRUE_(prime_f<1>( 43) == 193);
		TRUE_(prime_f<1>( 44) == 197);
		TRUE_(prime_f<1>( 45) == 199);
		TRUE_(prime_f<1>( 46) == 211);
		TRUE_(prime_f<1>( 47) == 223);
		TRUE_(prime_f<1>( 48) == 227);
		TRUE_(prime_f<1>( 49) == 229);
		TRUE_(prime_f<1>( 50) == 233);
		TRUE_(prime_f<1>( 51) == 239);
		TRUE_(prime_f<1>( 52) == 241);
		TRUE_(prime_f<1>( 53) == 251);
		TRUE_(prime_f<1>( 54) == 257);
		TRUE_(prime_f<1>( 55) == 263);
		TRUE_(prime_f<1>( 56) == 269);
		TRUE_(prime_f<1>( 57) == 271);
		TRUE_(prime_f<1>( 58) == 277);
		TRUE_(prime_f<1>( 59) == 281);
		TRUE_(prime_f<1>( 60) == 283);
		TRUE_(prime_f<1>( 61) == 293);
		TRUE_(prime_f<1>( 62) == 307);
		TRUE_(prime_f<1>( 63) == 311);
		TRUE_(prime_f<1>( 64) == 313);
		TRUE_(prime_f<1>( 65) == 317);
		TRUE_(prime_f<1>( 66) == 331);
		TRUE_(prime_f<1>( 67) == 337);
		TRUE_(prime_f<1>( 68) == 347);
		TRUE_(prime_f<1>( 69) == 349);
		TRUE_(prime_f<1>( 70) == 353);
		TRUE_(prime_f<1>( 71) == 359);
		TRUE_(prime_f<1>( 72) == 367);
		TRUE_(prime_f<1>( 73) == 373);
		TRUE_(prime_f<1>( 74) == 379);
		TRUE_(prime_f<1>( 75) == 383);
		TRUE_(prime_f<1>( 76) == 389);
		TRUE_(prime_f<1>( 77) == 397);
		TRUE_(prime_f<1>( 78) == 401);
		TRUE_(prime_f<1>( 79) == 409);
		TRUE_(prime_f<1>( 80) == 419);
		TRUE_(prime_f<1>( 81) == 421);
		TRUE_(prime_f<1>( 82) == 431);
		TRUE_(prime_f<1>( 83) == 433);
		TRUE_(prime_f<1>( 84) == 439);
		TRUE_(prime_f<1>( 85) == 443);
		TRUE_(prime_f<1>( 86) == 449);
		TRUE_(prime_f<1>( 87) == 457);
		TRUE_(prime_f<1>( 88) == 461);
		TRUE_(prime_f<1>( 89) == 463);
		TRUE_(prime_f<1>( 90) == 467);
		TRUE_(prime_f<1>( 91) == 479);
		TRUE_(prime_f<1>( 92) == 487);
		TRUE_(prime_f<1>( 93) == 491);
		TRUE_(prime_f<1>( 94) == 499);
		TRUE_(prime_f<1>( 95) == 503);
		TRUE_(prime_f<1>( 96) == 509);
		TRUE_(prime_f<1>( 97) == 521);
		TRUE_(prime_f<1>( 98) == 523);
		TRUE_(prime_f<1>( 99) == 541);
		TRUE_(prime_f<1>(100) == 547);
		TRUE_(prime_f<1>(101) == 557);
		TRUE_(prime_f<1>(102) == 563);
		TRUE_(prime_f<1>(103) == 569);
		TRUE_(prime_f<1>(104) == 571);
		TRUE_(prime_f<1>(105) == 577);
		TRUE_(prime_f<1>(106) == 587);
		TRUE_(prime_f<1>(107) == 593);
		TRUE_(prime_f<1>(108) == 599);
		TRUE_(prime_f<1>(109) == 601);
		TRUE_(prime_f<1>(110) == 607);
		TRUE_(prime_f<1>(111) == 613);
		TRUE_(prime_f<1>(112) == 617);
		TRUE_(prime_f<1>(113) == 619);
		TRUE_(prime_f<1>(114) == 631);
		TRUE_(prime_f<1>(115) == 641);
		TRUE_(prime_f<1>(116) == 643);
		TRUE_(prime_f<1>(117) == 647);
		TRUE_(prime_f<1>(118) == 653);
		TRUE_(prime_f<1>(119) == 659);
		TRUE_(prime_f<1>(120) == 661);
		TRUE_(prime_f<1>(121) == 673);
		TRUE_(prime_f<1>(122) == 677);
		TRUE_(prime_f<1>(123) == 683);
		TRUE_(prime_f<1>(124) == 691);
		TRUE_(prime_f<1>(125) == 701);
		TRUE_(prime_f<1>(126) == 709);
		TRUE_(prime_f<1>(127) == 719);
	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
