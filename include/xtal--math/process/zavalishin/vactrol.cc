#pragma once
#include "./any.cc"

#include "./prewarped.hh"
#include "./staged.hh"
#include "./gate.hh"

#include "./vactrol.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("vactrol")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;
	using U0_cue  = flow::cue_s<>;

	using U_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	U_alpha constexpr omega = 2*2*2*3*3*5*5;
	U_alpha constexpr   rho = 1;
	U_alpha constexpr    up = 1;
	U_alpha constexpr    dn = 0;

	/**/
	TRY_("vactrol: monophony")
	{
		using Z_filter = filter<>;
		using Z = any<Z_filter>;

		using refade_type = typename Z::refade_type;
		using redamp_type = typename Z::redamp_type;
		using reshape_type = typename Z::reshape_type;
		using   shape_type = typename Z::  shape_type;
		using   stage_type = typename Z::  stage_type;

		using Z_packet = flow::packet_t<stage_type, reshape_type>;
		using Z_event  = flow::cue_s<Z_packet>;
		using Z_cue    = flow::cue_s<>;

		using Z_process = confined_t<void
		,	prewarped<ordinal_constant_t<0>>, gate<1>
		,	typename redamp_type::template attend<>
		,	typename refade_type::template attend<>
		,	vactrol <>
		,	staged<-1>
		,	staged< 0>
		,	filter  <>
		>;
		using Z_processor = processor::monomer_t<Z_process
		//\
		,	U_slicer::template inqueue<stage_type, reshape_type>
		,	U_slicer::template inqueue<Z_packet>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= typename Z::   limit_type{0};
		z <<= typename Z::   order_type{2};
		z <<= typename Z::   patch_type{0};
		z <<= typename Z:: redamp_type{root_f<2>(2.F)};
		z <<= typename Z:: refade_type{0};

	//	z <<= reshape_type({0.25, one - 0.25});

		z <<= z_sample;
		z <<= z_resize;

		z >>= U_stage{-1};
		z >>= U0_cue{0x08}.then(Z_packet{ 0, shape_type{0.125, 0.25}});
		z >>= U0_cue{0x18}.then(Z_packet{-1, shape_type{0.500, 0.25}});

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(0 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08, 0x18);

		TRUE_(0 == z.efflux(z_cursor++));
//		TRUE_(1 == z.influx(occur::stage_f(-1)));
//
		echo_plot_<28>(z.store());

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
