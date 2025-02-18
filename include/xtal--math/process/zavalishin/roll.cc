#pragma once
#include "./any.cc"
#include "./prewarped.hh"
#include "./staged.hh"
#include "./trigger.hh"


#include "./roll.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("roll")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;
	using U0_cue  = flow::cue_s<>;

	using U_chunk = schedule::chunk_t<provision::spooled<extent_constant_t<0x10>>>;

	U_alpha constexpr omega = 2*3*5*5*7;
	U_alpha constexpr   rho = 1;
	U_alpha constexpr    up = 1;
	U_alpha constexpr    dn = 0;

	/**/
	TRY_("roll: monophony")
	{
		using U_value =             U_stage ;
		using U_event = flow::cue_s<U_value>;

		using  _process = any<filter<>>;
		using Z_process = confined_t<void
		,	prewarped<ordinal_constant_t<0>>
		,	trigger <>
		,	roll    <>
		,	staged<-1>
		,	staged< 0>
		,	typename _process::balance_type::template attend<>
		,	filter  <>
		>;
		using Z_processor = processor::monomer_t<Z_process
		,	U_chunk::template inqueue<U_value>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::bind_f(processor::let_f(f_));

		using recurve_type = typename _process::recurve_type;

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= typename _process::   limit_type{0};
		z <<= typename _process::   order_type{2};
		z <<= typename _process::   patch_type{0};
		z <<= typename _process:: damping_type{1};
		z <<= typename _process:: balance_type{0.5};

		z <<= recurve_type({+0.25, one - 0.25});

		z <<= z_sample;
		z <<= z_resize;
		z <<= U_stage(     -1);

		z <<= U_event(0x08,  0);
		z <<= U_event(0x10,  0);
		z <<= U_event(0x28, -1);
	//	z <<= U_event(0x38,  0);

		echo_rule_<25>('=');
		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.efflux(occur::stage_f(-1)));

		echo_plot_<25>(z.store());

		echo_rule_<25>();
		TRUE_(0 == z.efflux(z_cursor++));
//		TRUE_(1 == z.efflux(occur::stage_f(-1)));
//
		echo_plot_<25>(z.store());
//		echo_rule_<25>();

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
