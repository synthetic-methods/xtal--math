#pragma once
#include "./any.cc"
#include "./intake.hh"
#include "./retake.hh"
#include "../../provision/prewarping.hh"


#include "./vactrol.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("vactrol")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using W_alpha = atom::math::dot_t<U_alpha[2]>;
	using Z_slice = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("vactrol: monophony")
	{
		//\
		using E_def = filter<U_alpha[2], union ENV>;
		using E_def = filter<>;
		using E_etc = occur::context_t<E_def>;
		using E_pkt = flow::packet_t<typename E_etc::stage_type, typename E_etc::reshape_type>;
		using E_prx = confined_t<void
		,	provision::math::prewarping< 0>
		,	intake< 1>
		,	retake< 0>
		,	retake<-1>
		,	typename E_etc::redamp_type::template   attend<>
		,	typename E_etc::refade_type::template   attend<>
		,	typename E_etc::             template   attach<>
		,	typename E_etc::             template dispatch<>
		,	vactrol<>
		,	E_def
		,	provision::math::saturation<identity>
		>;
		using E_pxr = processor::monomer_t<E_prx
		//\
		,	Z_slice::template accept<stage_type, typename E_etc::reshape_type>
		,	Z_slice::template accept<E_pkt>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr e_omega = 2*2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = E_pxr::bind_f(processor::let_f(e_omega));
		z <<= typename E_etc::  order_type{2};
		z <<= typename E_etc:: redamp_type{root_f<2>(2.F)};
		z <<= typename E_etc:: refade_type{0};
	//	z <<= typename E_etc::reshape_type({0.25, one - 0.25});

		z <<= z_sample;
		z <<= z_resize;
		z >>= typename E_etc::stage_type{-1};

		z >>= flow::cue_f(0x08).then(E_pkt{ 0, typename E_etc::shape_type{0.125, 0.25}});
		z >>= flow::cue_f(0x18).then(E_pkt{-1, typename E_etc::shape_type{0.500, 0.25}});

		echo_("\nvactrol: monophony");
	//	echo_rule_<28>("\u2500");

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
