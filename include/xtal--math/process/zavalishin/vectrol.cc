#pragma once
#include "./any.cc"

#include "./prewarped.hh"
#include "./staged.hh"
#include "./gate.hh"

#include "./vectrol.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("vectrol")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;
	using U0_cue  = flow::cue_s<>;

	using U_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	U_alpha constexpr omega = 2*2*3*3*5*5;
	U_alpha constexpr   rho = 1;
	U_alpha constexpr    up = 1;
	U_alpha constexpr    dn = 0;

	/**/
	TRY_("vectrol: monophony")
	{
		using U_value =             U_stage ;
		using U_event = flow::cue_s<U_value>;

		using A_filter = filter<>;

		using A = any<A_filter>;
		using Z_process = confined_t<void
		,	prewarped<ordinal_constant_t<0>>, gate<0>
		,	typename A::damping_type::template attend<>
		,	typename A::balance_type::template attend<>
		,	vectrol    <>
		,	staged<-1>
		,	staged< 0>
		,	A_filter
		>;
		using Z_processor = processor::monomer_t<Z_process
		,	U_slicer::template inqueue<U_value>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::bind_f(processor::let_f(f_));

		using recurve_type = typename A::recurve_type;

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= typename A::   limit_type{0};
		z <<= typename A::   order_type{2};
		z <<= typename A::   patch_type{0};
		z <<= typename A:: damping_type{1};
		z <<= typename A:: balance_type{0.5};

		z <<= recurve_type({0.25, one - 0.25});

		z <<= z_sample;
		z <<= z_resize;
		z <<= U_stage(     -1);

		z <<= U_event(0x08,  0);
		z <<= U_event(0x18,  0);
		z <<= U_event(0x28, -1);
	//	z <<= U_event(0x38,  0);

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.efflux(occur::stage_f(-1)));

		echo_rule_<25>('=');
		echo_plot_<25>(z.store());

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(1 == z.efflux(occur::stage_f(-1)));

		echo_rule_<25>('-');
		echo_plot_<25>(z.store());

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
