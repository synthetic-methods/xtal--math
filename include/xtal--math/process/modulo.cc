#pragma once
#include "./any.cc"





#include "./modulo.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("math")
{
	using U_fit   = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;
	using U_phi   = atom::math::phason_t<U_alpha[2]>;

	using W_alpha = atom::couple_t<U_alpha[2]>;
	using W_aphex = atom::couple_t<U_aphex[2]>;

	TRY_("modulo > 0")
	{
		TRUE_(modulo_f<3U>(0) == 0);
		TRUE_(modulo_f<3U>(4) == 1);

	}
	TRY_("modulo < 0")
	{
		U_phi q0{0.45, 0.45};
		U_phi q1(modulo_f<-1>(q0));
		TRUE_(check_f<-48>(0.05, modulo_f<-1>(0.45)));
		TRUE_(check_f<-48>(0.05, q1(0)));
		TRUE_(check_f<-18>(0.45, q1(1)));

		TRUE_(modulo_f<0b111U>(10) == 10%7);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
