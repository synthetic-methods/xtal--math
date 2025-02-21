#pragma once
#include "./any.cc"

#include "./prewarped.hh"



#include "./staged.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("staged")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;
	using U0_cue   = flow::cue_s<>;

	using U_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	U_alpha constexpr omega = 2*2*2*3*5*5*7;
	U_alpha constexpr   rho = 1;
	U_alpha constexpr    up = 1;
	U_alpha constexpr    dn = 0;

	/*/
	TRY_("staged: ?")
	{
		TRUE_(true);

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
