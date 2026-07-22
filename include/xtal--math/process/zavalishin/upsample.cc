#pragma once
#include "./any.cc"

#include "./reuse.hh"



#include "./upsample.hh"
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("upsample")
{
	using U_fit   = bond::fit<>;
	using U_sigma = U_fit::sigma_type;
	using U_delta = U_fit::delta_type;
	using U_alpha = U_fit::alpha_type;
	using U_aphex = U_fit::aphex_type;

	using W_alpha = atom::math::dot_t<U_alpha[2]>;
	using Z_slice = schedule::slicer_t<scheme::spooled<extent_constant_t<0x10>>>;

	using X_resample = occur::resample_t<>;
	using A_resample = typename X_resample::template attach<>;
	//\
	using U_shape = atom::                   couple_t<U_alpha[2]>;
	using U_shape = atom:: quantity_plus_multiplies_t<U_alpha[2]>;
	using U_coeff = atom::math::                dot_t<U_alpha[2]>;

	using X_shape = occur::reinferred_t<U_shape, union SHAPE>;
	using X_coeff = occur::reinferred_t<U_coeff, union COEFF>;

	using X_stage = occur::stage_t<>;
	using Y_trig  = pulse_t< 0>;
	using Y_gate  = pulse_t< 1>;
	using Y_hold  = pulse_t<-1>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/*/
	TRY_("upsample: check")
	{
		TRUE_(true);

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
