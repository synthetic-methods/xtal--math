#pragma once
#include "./any.cc"
#include "./gate.hh"
#include "./staged.hh"
#include "./prewarped.hh"


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
		using E_env = any_t<E_def>;
		using E_eve = flow::packet_t<typename E_env::stage_type, typename E_env::reshape_type>;
		using E_pro = confined_t<void
		,	prewarped<_0>, gate<1>
		,	staged<-1>
		,	staged< 0>
		,	typename E_env::redamp_type::template   attend<>
		,	typename E_env::refade_type::template   attend<>
		,	typename E_env::             template   attach<>
		,	typename E_env::             template dispatch<>
		,	vactrol<>
		,	E_def
		>;
		using E_prx = processor::monomer_t<E_pro
		//\
		,	Z_slice::template inqueue<stage_type, typename E_env::reshape_type>
		,	Z_slice::template inqueue<E_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr e_omega = 2*2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = E_prx::bind_f(processor::let_f(e_omega));
		z <<= typename E_env::  limit_type{0};
		z <<= typename E_env::  order_type{2};
		z <<= typename E_env::  patch_type{0};
		z <<= typename E_env:: redamp_type{root_f<2>(2.F)};
		z <<= typename E_env:: refade_type{0};
	//	z <<= typename E_env::reshape_type({0.25, one - 0.25});

		z <<= z_sample;
		z <<= z_resize;
		z >>= typename E_env::stage_type{-1};

		z >>= flow::cue_f(0x08).then(E_eve{ 0, typename E_env::shape_type{0.125, 0.25}});
		z >>= flow::cue_f(0x18).then(E_eve{-1, typename E_env::shape_type{0.500, 0.25}});

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
