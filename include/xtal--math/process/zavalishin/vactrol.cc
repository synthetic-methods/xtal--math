#pragma once
#include "./any.cc"

#include "./reuse.hh"



#include "./vactrol.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("vactrol")
{
	using U_fit   = bond::fit<>;
	using U_sigma = U_fit::sigma_type;
	using U_delta = U_fit::delta_type;
	using U_alpha = U_fit::alpha_type;
	using U_aphex = U_fit::aphex_type;

	using W_alpha = atom::math::dot_t<U_alpha[2]>;
	using Z_slice = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	using X_resample = occur::resample_t<>;
	using A_resample = typename X_resample::template attach<>;

	using U_shape = atom::      couple_t<U_alpha[2]>;
	using U_coeff = atom::math::   dot_t<U_alpha[2]>;

	using X_shape = occur::reinferred_t<U_shape, union SHAPE>;
	using X_coeff = occur::reinferred_t<U_coeff, union COEFF>;

	using X_stage = occur::stage_t<>;
	using Y_trig  = pulse_t< 0>;
	using Y_gate  = pulse_t< 1>;
	using Y_hold  = pulse_t<-1>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("vactrol: monophonic gate")
	{
		//\
		using F_env = filter<>;
		using F_env = filter<U_alpha[2], union ENV>;
		using L_env = occur::codex_t<F_env>;
		using X_env = flow::packet_t<X_stage, X_shape>;
		using Y_env = confined_t<void
		,	reuse< 0, -1>
		,	coefficient_t<X_coeff> ::   attach <>
		,	Y_gate                 ::   infix  <>
		,	X_shape                ::   affix  <>
		,	L_env                  ::   attach <>
		,	L_env                  :: dispatch <>
		,	vactrol<1>
		,	F_env
		,	provision::math::zavalishin::shaped<identity>
		>;
		using Z_env = processor::monomer_t<Y_env
		,	Z_slice::template suspend<X_env>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		//\
		using Y_phi = phasor_t<U_alpha[2], typename X_resample::template attach<>>;
		using Y_phi = confined_t<phasor<U_alpha[2]>, typename X_resample::attach<>>;
		using Z_phi = processor::monomer_t<Y_phi>;

		U_alpha constexpr e_omega = 2*2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = Z_env::bind_f(Z_phi::bind_f(processor::let_f(e_omega)));
		z <<= typename L_env::order_attribute{2};
		z <<= X_coeff{1, 0};
		z >>= X_stage{  -1};
		z <<= z_sample;
		z <<= z_resize;

		z >>= flow::cue_f(0x08).then(X_env{ X_stage{0}, X_shape{-0.50, -1.00}});
		z >>= flow::cue_f(0x18).then(X_env{ X_stage{1}, X_shape{ 1.00, -0.00}});

		echo_("\nvactrol: monophony");

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(0 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08, 0x18);

		TRUE_(0 == z.efflux(z_cursor++));
//		TRUE_(1 == z.influx(occur::stage_f(-1)));
//
		echo_plot_<28>(z.store());

	}
	/***/
	/**/
	TRY_("vactrol: monophonic trigger")
	{
		//\
		using F_env = filter<>;
		using F_env = filter<U_alpha[2], union ENV>;
		using L_env = occur::codex_t<F_env>;
		using X_env = flow::packet_t<X_stage, X_shape>;
		using Y_env = confined_t<void
		,	reuse< 0, -1>
		,	coefficient_t<X_coeff> ::   attach <>
		,	Y_trig                 ::   infix  <>
		,	X_shape                ::   affix  <>
		,	L_env                  ::   attach <>
		,	L_env                  :: dispatch <>
		,	vactrol<0>
		,	F_env
		,	provision::math::zavalishin::shaped<identity>
		>;
		using Z_env = processor::monomer_t<Y_env
		,	Z_slice::template suspend<X_env>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		//\
		using Y_phi = phasor_t<U_alpha[2], typename X_resample::template attach<>>;
		using Y_phi = confined_t<phasor<U_alpha[2]>, typename X_resample::attach<>>;
		using Z_phi = processor::monomer_t<Y_phi>;

		U_alpha constexpr e_omega = 1*2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = Z_env::bind_f(Z_phi::bind_f(processor::let_f(e_omega)));
		z <<= typename L_env::order_attribute{2};
		z <<= X_coeff{1, 0};
		z >>= X_stage{  -1};
		z <<= z_sample;
		z <<= z_resize;

		z >>= flow::cue_f(0x08).then(X_env{ X_stage{0}, X_shape{0.125,  1.00}});
	//	z >>= flow::cue_f(0x18).then(X_env{ X_stage{1}, X_shape{1.000, -0.00}});

		echo_("\nvactrol: monophony");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(0 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08);

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
