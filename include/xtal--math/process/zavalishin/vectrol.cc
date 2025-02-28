#pragma once
#include "./any.cc"
#include "./gate.hh"
#include "./staged.hh"
#include "./prewarped.hh"


#include "./vectrol.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("vectrol")
{
	using T_alpha   = typename bond::fit<>::alpha_type;
	using T_aphex   = typename bond::fit<>::aphex_type;
	//\
	using A_filter  = filter<T_alpha[2], union RING>;
	using A_filter  = filter<>;
	using W_filter  = any_t<A_filter>;
	using U_limit   = typename W_filter::   limit_type;
	using U_order   = typename W_filter::   order_type;
	using U_patch   = typename W_filter::   patch_type;
	using U_stage   = typename W_filter::   stage_type;
	using U_shape   = typename W_filter::   shape_type;
	using U_refade  = typename W_filter::  refade_type;
	using U_redamp  = typename W_filter::  redamp_type;
	using U_reshape = typename W_filter:: reshape_type;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	using Z_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	/**/
	TRY_("vectrol: monophony")
	{
		using U_event =             U_stage ;
		using Z_event = flow::cue_s<U_event>;

		using Z_process = confined_t<void
		,	prewarped<_0>, gate<0>
		,	staged<-1>
		,	staged< 0>
		,	typename U_redamp::template   attend<>
		,	typename U_refade::template   attend<>
		,	typename W_filter::template   attach<>
		,	typename W_filter::template dispatch<>
		,	vectrol<>
		,	A_filter
		>;
		using Z_processor = processor::monomer_t<Z_process
		,	Z_slicer::template inqueue<U_event>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		T_alpha constexpr omega = 2*2*3*3*5*5;
		auto z = Z_processor::bind_f(processor::let_f(omega));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= U_limit  {0};
		z <<= U_order  {2};
		z <<= U_patch  {0};
		z <<= U_redamp {1};
		z <<= U_refade {0.5};

		z <<= U_reshape{U_shape{0.125, one - 0.125}};

		z <<= z_sample;
		z <<= z_resize;

		z >>= U_stage{-1};
		z >>= Z_event(0x08,  0);
		z >>= Z_event(0x18,  0);
		z >>= Z_event(0x28, -1);
	//	z >>= Z_event(0x38,  0);

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08, 0x18);

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(1 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08);

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
