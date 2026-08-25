#pragma once
#include "./any.cc"
#include "../../occur/all.hh"
#include "../../scheme/all.hh"
#include "./reuse.hh"
#include "./vactrol.hh"

#include "./filter.hh"
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test::XTAL_NUM
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
using namespace xtal::process::math::zavalishin;

template <class U>
using arguments_x = atom::bucket_t<U[2]>;

struct filter_parameters
{
	template <class S>
	class subtype : public S
	{
	public:
		using S::S;
		using test_type = bool;

	};
};


////////////////////////////////////////////////////////////////////////////////

TAG_("filter")
{
	using U_fit   = typename bond::fit<>;
	using U_alpha = U_fit::alpha_type;
	using W_alpha = atom::math::dot_t<U_alpha[2]>;
	using Z_slice = schedule::slicer_t<scheme::spooled<extent_constant_t<0x10>>>;

	using U_quartz = occur::quartz_t<>;
	using U_rewind = occur::rewind_t<occur::stage_t<>>;


	using U_coeff  = atom::math::dot_t<U_alpha[2]>;
	using X_coeff  = occur::inferred_t<U_coeff, union COEFF>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("filter: 1D instantiation")
	{
		using R_def = filter<>;
		using R_etc = process::occurrence_t<R_def>;
		using R_prx = confined_t<void
		,	per_t<U_quartz>            :: transfix<1>
		,	U_rewind                   ::   attach< >
		,	coefficient_t<X_coeff>     ::   attach< >
	//	,	reuse< 0>
		,	R_etc                      ::   attach< >
		,	R_etc                      :: dispatch< >
		,	R_def
		,	scheme::math::zavalishin::distorted<identity>
		>;
		using R_pxr = processor::monomer_t<R_prx>;

		R_prx svf{};

		auto constexpr N0_sample_rate  =  (2*2)*(3*3)*(5*5)*(7*7);    // 44100
		auto constexpr N0_filter_rate  =  (2*2)*(3*3)*(5*5)*(7);      //  6300
		auto constexpr N1_filter_rate  =  U_fit::ratio_f(N0_filter_rate, N0_sample_rate);
		auto constexpr K1_filter_rate  =  pade::tangy_f<1>(                       N1_filter_rate);
	//	auto constexpr K1_growth_rate  =  pade::tangy_f<1>(U_fit::ratio_f(1, 4) - N1_filter_rate);
		auto constexpr K1_growth_rate  =  (one - K1_filter_rate)/(one + K1_filter_rate);

		svf <<= occur::quartz_f(N0_sample_rate);
		svf <<= typename R_etc::order_attribute{1};
		svf <<= X_coeff{1, 0};

		U_alpha constexpr omega = N0_filter_rate;
		U_alpha constexpr   rho = 1;
		U_alpha constexpr    up = 1;
		U_alpha constexpr    dn = 0;

		U_alpha v1{one};
		U_alpha v0{one};

		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);// TRUE_(check_f<- 0>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));

		v1 = one - v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);

	}
	TRY_("filter: 2D instantiation")
	{
		using R_def = filter<>;
		using R_etc = process::occurrence_t<R_def>;
		using R_prx = confined_t<void
		,	per_t<U_quartz>            :: transfix<1>
		,	U_rewind                   ::   attach< >
		,	coefficient_t<X_coeff>     ::   attach< >
	//	,	reuse< 0>
		,	R_etc                      ::   attach< >
		,	R_etc                      :: dispatch< >
		,	R_def
		,	scheme::math::zavalishin::distorted<identity>
		>;
		using R_pxr = processor::monomer_t<R_prx>;

		R_prx svf{};

		auto constexpr N0_sample_rate  =  (2*2)*(3*3)*(5*5)*(7*7);    // 44100
		auto constexpr N0_filter_rate  =  (2*2)*(3*3)*(5*5)*(7);      //  6300
		auto constexpr N1_filter_rate  =  U_fit::ratio_f(N0_filter_rate, N0_sample_rate);
		auto constexpr K1_filter_rate  =  pade::tangy_f<1>(                       N1_filter_rate);
	//	auto constexpr K1_growth_rate  =  pade::tangy_f<1>(U_fit::ratio_f(1, 4) - N1_filter_rate);
		auto constexpr K1_growth_rate  =  (one - K1_filter_rate)/(one + K1_filter_rate);

		svf <<= occur::quartz_f(N0_sample_rate);
		svf <<= typename R_etc::order_attribute{2};
		svf <<= X_coeff{1, 1};

		U_alpha constexpr omega = N0_filter_rate;
		U_alpha constexpr   rho = 1;
		U_alpha constexpr    up = 1;
		U_alpha constexpr    dn = 0;

		U_alpha v1{one};
		U_alpha v0{one};

		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);// TRUE_(check_f<- 0>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));
		v1 =       v0; v0 = one - svf(up, omega, rho); TRUE_(v0 < v1);   TRUE_(check_f<-22>(K1_growth_rate, v0/v1));

		v1 = one - v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);
		v1 =       v0; v0 =       svf(dn, omega, rho); TRUE_(v1 > v0);

	}
	TRY_("filter parameterization")
	{
		using Y_ramp = confined_t<void
		,	typename per_t<U_quartz>::template transfix<1>
		//\
		,	filter<union RAMP>
		,	filter<U_alpha[2], union RAMP>
		>;
		using Y_ring = confined_t<void
		,	typename per_t<U_quartz>::template transfix<1>
		//\
		,	filter<union RING>
		,	filter<U_alpha[2], union RING>
		>;
		static_assert(different_q<typename Y_ramp::order_attribute, typename Y_ring::order_attribute>);

	};
	/***/
}

////////////////////////////////////////////////////////////////////////////////

TAG_("filter-ring")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using T_aphex = typename bond::fit<>::aphex_type;
	using Z_slice = schedule::slicer_t<scheme::spooled<extent_constant_t<0x10>>>;

	using U_quartz = occur::quartz_t<>;
	using U_rewind = occur::rewind_t<>;
	using U_stage  = occur::stage_t<>;

	using U_coeff  = atom::math::dot_t<U_alpha[2]>;
	using X_coeff  = occur::  inferred_t<U_coeff, union COEFF>;
	using X_gain   = occur::reinferred_t<U_alpha, union GAIN >;

	using Y_trig  = pulse_t< 0>;
	using Y_gate  = pulse_t< 1>;
	using Y_hold  = pulse_t<-1>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("filter-ring monophony")
	{
		using R_def = filter<U_alpha[2], union RING>;
		using R_etc = process::occurrence_t<R_def>;
		using R_eve = flow::packet_t<U_stage, typename R_etc::damp_parameter>;
		using R_prx = confined_t<void
		,	per_t<U_quartz>::transfix<0>
		//\
		,	reuse<0, -1>
		,	reuse<   -1>
		,	U_rewind               ::   attach< >
		,	coefficient_t<X_coeff> ::   attach< >
	//	,	Y_hold                 ::   prefix< >
		,	X_gain                 ::   prefix< >
		,	  pulse_t<-1>          :: transfix<0>
		,	vactrol_t<-1>          :: transfix<0>
		,	R_etc::damp_parameter  ::   suffix< >
		,	R_etc::                     attach< >
		,	R_etc::                   dispatch< >
		,	R_def
		>;
		using R_pxr = processor::monomer_t<R_prx
		,	Z_slice::template suspend<R_eve>
		,	scheme::stored  <null_type[0x100]>
		,	scheme::spooled <null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 2*2*2*3*5*5*7;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::quartz_f  (44100);

		auto z = R_pxr::bind_f(processor::let_f(r_omega));
		z <<= typename R_etc::order_attribute{2};
		z <<= typename R_etc:: damp_parameter{1};
		z <<= X_coeff{0, 1};
		z <<= X_gain {1};

		z <<= typename occur::rewind_t<>{0};

		z <<= z_sample;
		z <<= z_resize;

		z >>= occur::stage_f(-1);
		z >>= flow::cue_f(0x08).then(R_eve{ 0, 0.0F});
		z >>= flow::cue_f(0x10).then(R_eve{ 0, 0.0F});
		//\
		z >>= flow::cue_f(0x20).then(R_eve{ 1, 0.5F});
		z >>= flow::cue_f(0x21).then(R_eve{ 1, 0.5F});

		echo_("\nfilter-ring: monophony");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(1 == z.influx(occur::stage_f( 0)));// Would be unchanged...
		TRUE_(0 == z.influx(occur::stage_f( 1)));
		TRUE_(0 == z.influx(occur::stage_f(-1)));
		echo_plot_<28>(z.store(), 0x08, 0x10);

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.influx(occur::stage_f( 0)));
		TRUE_(1 == z.influx(occur::stage_f( 1)));// Would be unchanged...
		TRUE_(0 == z.influx(occur::stage_f(-1)));
		//\
		echo_plot_<28>(z.store(), 0x00);
		echo_plot_<28>(z.store(), 0x01);

	}
	/***/
	/**/
	TRY_("filter-ring polyphony")
	{
		using R_def = filter<U_alpha[2], union RING>;
		using R_etc = process::occurrence_t<R_def>;
		using R_eve = flow::key_s<U_stage>;

		using R_prx = confined_t<void
	//	,	reuse<0, -1>
		,	reuse<   -1>
		,	per_t<U_quartz>        :: transfix<0>
		,	coefficient_t<X_coeff> ::   attach< >
		,	U_rewind               ::   attach< >
	//	,	Y_gate                 ::   prefix< >
		,	X_gain                 ::   prefix< >
		,	  pulse_t< 1>          :: transfix< >
		,	vactrol_t< 1>          :: transfix< >
		,	U_stage::template assignment<typename R_etc::damp_parameter>
		,	R_etc::damp_parameter  ::   suffix< >
		,	R_etc                  ::   attach< >
		,	R_etc                  :: dispatch< >
		,	R_def
		>;
		using R_pxy = processor::polymer_t<R_prx
		,	Z_slice::template suspend<R_eve>
		,	scheme::stored <null_type[0x100]>
		,	scheme::spooled<null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 59*61;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::quartz_f  (44100);

		auto z = R_pxy::bind_f(processor::let_f(r_omega));
		z <<= typename R_etc::order_attribute{2};
		z <<= typename R_etc:: damp_parameter{1};
		z <<= X_coeff{0, 1};
		z <<= X_gain{0.707};
		z <<= typename occur::rewind_t<>{0};

		z <<= flow::assign_f(U_stage{ 0}) << typename R_etc::damp_parameter{0.0000F};
		z <<= flow::assign_f(U_stage{ 1}) << typename R_etc::damp_parameter{0.0625F};
		z <<= flow::assign_f(U_stage{-1}) << typename R_etc::damp_parameter{0.5000F};

		z <<= z_sample;
		z <<= z_resize;

		z.lead() >>= U_stage{-1};

		z >>= flow::cue_f(0x04).then(R_eve{69, 0});
		z >>= flow::cue_f(0x10).then(R_eve{69, 1});
		z >>= flow::cue_f(0x18).then(R_eve{69, 0});
		z >>= flow::cue_f(0x28).then(R_eve{69, 1});
	//	z >>= flow::cue_f(0x38).then(R_eve{69,-1});
	//	z >>= flow::cue_f(0x40).then(R_eve{69, 0});
	//	z >>= flow::cue_f(0x50).then(R_eve{69, 1});

		TRUE_(0 == z.ensemble().size());
		echo_("\nfilter-ring: polyphony");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(2 == z.ensemble().size());// Still decaying!
		echo_plot_<28>(z.store(), 0x04, 0x10, 0x18);

	//	z >>= flow::cue_f(0x08).then(R_eve{69, 1});// Inlined below...
		z >>= flow::cue_f(0x20).then(R_eve{69,-1});// Inlined from above...

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(2 >= z.ensemble().size());// Still decaying?
		echo_plot_<28>(z.store(), 0x08);

	//	z >>= flow::cue_f(0x00).then(R_eve{69, 0});// Outlined above...
	//	z >>= flow::cue_f(0x08).then(R_eve{69,-1});// Outlined above...

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(1 == z.ensemble().size());// Still decaying...
		echo_plot_<28>(z.store(), 0x00);

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.ensemble().size());
	//	echo_plot_<28>(z.store());

	}
	/***/
}


////////////////////////////////////////////////////////////////////////////////

TAG_("vectrol")
{
	using U_alpha  = typename bond::fit<>::alpha_type;
	using W_alpha  = atom::math::dot_t<U_alpha[2]>;
	using Z_slice  = schedule::slicer_t<scheme::spooled<extent_constant_t<0x10>>>;

	using U_quartz = occur::quartz_t<>;
	using U_stage  = occur::stage_t<>;

	using U_coeff  = atom::math::dot_t<U_alpha[2]>;
	using X_coeff  = occur::  inferred_t<U_coeff, union COEFF>;
	using X_gain   = occur::reinferred_t<U_alpha, union GAIN >;
	using X_amph   = occur::reinferred_t<U_alpha, union AMPH >;

	using Y_trig   = pulse_t< 0>;
	using Y_gate   = pulse_t< 1>;
	using Y_hold   = pulse_t<-1>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("vectrol: monophony")
	{
		//\
		using S_content = filter<>;
		using S_content = filter<U_alpha[2], union ENV>;
		//\
		using S_meta = confined_t<S_content>;
		using S_meta = process::occurrence_t<S_content>;

		using S_damp_   = occur::math::zavalishin::probe_t<typename S_meta::codata_type>;
	//	using S_damp    = typename S_meta::damp_parameter;

		using S_order   = typename S_meta::order_attribute;

		using S_process = confined_t<void
		,	reuse<0, -1>
		,	coefficient_t<X_coeff> ::   attach< >
		,	per_t<U_quartz>        :: transfix<0>
	//	,	Y_trig                 ::   prefix< >
		,	X_gain                 ::   prefix< >
		,	  pulse_t< 0>          :: transfix< >
		,	vactrol_t< 0>          :: transfix< >
		,	S_damp_                ::   suffix< >
	//	,	S_damp                 ::   suffix< >
		,	S_meta                 ::   attach< >
		,	S_meta                 :: dispatch< >
		,	S_content
		>;
		using S_processor = processor::monomer_t<S_process
		,	Z_slice::template suspend<U_stage>
		,	scheme::stored  <null_type[0x100]>
		,	scheme::spooled <null_type[0x100]>
		>;

		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::quartz_f  (44100);

		auto z = S_processor::bind_f(processor::let_f(e_omega));
		z <<= S_order{2};
		z <<= X_gain{1};
		z <<= S_damp_{one - 0.0, zero - 1.0};
		z <<= X_coeff{1, 1};

		z <<= z_sample;
		z <<= z_resize;
		z >>= U_stage{-1};

		z >>= flow::cue_f(0x08).then(U_stage{ 0});
		z >>= flow::cue_f(0x18).then(U_stage{ 0});
		z >>= flow::cue_f(0x28).then(U_stage{-1});
	//	z >>= flow::cue_f(0x38).then(U_stage{ 0});

		echo_("\nvectrol: monophony");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08, 0x18);

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(1 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08);

	}
	/***/
	/**/
	TRY_("vectrol: patch")
	{
		using S_content = filter<U_alpha[2], union ENV>;
		using S_meta    = process::occurrence_t<S_content>;
		using S_damp_   = occur::math::zavalishin::probe_t<typename S_meta::codata_type>;
		using S_damp    = typename S_meta::  damp_parameter;
		using S_order   = typename S_meta:: order_attribute;

		using S_process = confined_t<void
		,	reuse< 0, -1>
	//	,	process::lift<W_alpha>
		,	typename per_t<U_quartz> :: transfix<0>
	//	,	Y_trig                   ::   prefix< >
		,	X_gain                   ::   prefix< >
		,	  pulse_t< 0>            :: transfix< >
		,	vactrol_t< 0>            :: transfix< >
		,	typename S_damp_         ::   suffix< >
		,	typename S_meta          ::   attach< >
		,	typename S_meta          :: dispatch< >
		,	S_content
		>;
		//\
		using S_processor = processor::conferred_t<S_process
		using S_processor = processor::monomer_t<S_process
	//	,	Z_slice::template suspend<U_stage>
		,	scheme::stored  <unit_type[0x100]>
		,	scheme::spooled <null_type[0x100]>
		>;

		using T_content =  filter<U_alpha[2], union RING>;
		using T_meta    =  process::occurrence_t<T_content>;
	//	using T_damp_   =  occur::math::zavalishin::probe_t<typename T_meta::codata_type>;
		using T_damp    =  typename T_meta:: damp_parameter;

		using T_order   =  typename T_meta::order_attribute;
		using T_event   =  flow::packet_t<U_stage, T_damp>;
		using T_dummy   =  occur::inferred_t<U_alpha, union DUMMY>;

		//\
		using X_vector  =  atom::bucket_t<U_alpha, W_alpha>;
	//	using X_vector  =  flow::packed_t<U_alpha, W_alpha, union SHLEEM>;
	//	using X_matrix  =  atom::bucket_t<X_vector[2]>;

		using X_row     =  atom::bucket_t<U_alpha, W_alpha>;
		using X_ray     =  atom::bucket_t<null_type[2]>;
		using X_matrix  =  bond::compose_s<X_row, atom::bubble<X_ray::template recast_t>>;
		using X_patch   =  patch_t<X_matrix>;

	//	using Vec_one     =  occur::reinferred_t<U_alpha, union ONE>;
		using Var_pulse   =  occur::reinferred_t<float, union PULSE>;
		using Var_pitch   =  occur::reinferred_t<float, union PITCH>;
		using Var_head    =  occur::reinferred_t<float, union HEAD>;
		using Var_tail    =  occur::reinferred_t<float, union TAIL>;
		using Var_dyn     =  float;

		using Vec_dyn     =  flow::packed_t<Var_pulse, Var_pitch, Var_head, Var_tail, union DYN>;
		using Vtx_dyn     =  atom::bucket_t<Var_dyn[4], union DYN>;
		//\
		using Vxx_dyn     =                 Vtx_dyn;
		using Vxx_dyn     =                 Vec_dyn;
		using Mxx_dyn     =  atom::bucket_t<Vxx_dyn[4]>;
		using Pxx_dyn     =         patch_t<Mxx_dyn>;

		using Evt_dyn     =  atom::bucket_t<float[4], union DYN>;// one, key, vel, vax
	//	TODO: Fix autoconversion issue when underlying is mismatched.

		using X_process = confined_t<void
		,	Pxx_dyn                   ::   deflow<Evt_dyn>
	//	,	Pxx_dyn                   ::   reflow<Vec_dyn>
		,	Pxx_dyn                   ::   reflow< >
		,	Pxx_dyn                   ::   attach< >
		,	X_patch                   ::   rewire< >
		,	X_patch                   ::   attach< >
		,	reuse< 0>
		,	coefficient_t<X_coeff>    ::   attach< >
		,	coefficient_t<>           ::    unfix< >
		,	per_t<U_quartz>           :: transfix<0>
	//	,	Y_gate                    ::   prefix< >
		,	X_amph                    ::   prefix< >
		,	  pulse_t< 0>             :: transfix< >
		,	vactrol_t< 0>             :: transfix< >
	//	,	T_damp_                   ::   suffix< >
		,	T_damp                    ::   suffix< >
		,	T_meta                    ::   attach< >
		,	T_meta                    :: dispatch< >
		,	T_content
		>;
		using X_processor = processor::monomer_t<X_process
		,	Z_slice::template suspend<T_event>
		,	scheme::stored  <null_type[0x100]>
	//	,	scheme::spooled <null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 3*3*3*5*5*5;
		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::quartz_f  (44100);

		auto _1 = processor::let_f(U_alpha{one});
		auto _e = processor::let_f(e_omega);
		auto _r = processor::let_f(r_omega);
		auto _x = S_processor::bind_f(_e);
		//\
		auto _y = X_processor::bind_f(_1, _x);
		auto _y = X_processor::bind_f(_1, S_processor::bind_f(_e));
		//\
		auto _w = occur::math::dent_s<Mxx_dyn, 0>({0x0, 0x1, 0x2, 0x3});// Vtx
		auto _w = occur::math::dent_s<Mxx_dyn, 0> {0x0, 0x1, 0x2, 0x3} ;// Vec
		auto _v = occur::math::dent_s<Vec_dyn   > {0x0, 0x1, 0x2, 0x3};

		_y <<= occur::math::dent_s<Mxx_dyn , 0>{ 0.f,  1.f,  2.f,  3.f};
		_y <<= occur::math::dent_s<Mxx_dyn , 1>{ 4.f,  5.f,  6.f,  7.f};
		_y <<= occur::math::dent_s<Mxx_dyn , 2>{ 8.f,  9.f, 10.f, 11.f};
		_y <<= occur::math::dent_s<Mxx_dyn , 3>{12.f, 13.f, 14.f, 15.f};
		_y <<= Evt_dyn{0.88, 0.66, 0.44, 0.22};

		_y <<= occur::math::dent_s<X_matrix, 1>({{1.00, 1.00}, {1111, 1111}});
		_y <<= occur::math::dent_s<X_matrix, 0> { 0.00       ,  r_omega    } ;

		_y <<= S_order {2};
		_y <<= X_gain  (1);
		_y <<= X_amph  (1);
		_y <<= S_damp_ {one - 0.0, zero - 1.0};
	//	_y <<= S_damp  {1.0};

		_y <<= T_order {2};
	//	_y <<= T_damp_ {1};
		_y <<= T_damp  {1};
		_y <<= X_coeff {0, 1};

		_y <<= z_sample;
		_y <<= z_resize;
		_y >>= U_stage{-1};

		_y <<= z_resize;

		_y >>= flow::cue_f(0x08).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x08).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 1});
	//	_y >>= flow::cue_f(0x18).then(U_stage{ 1});

		_y >>= flow::cue_f(0x18).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x18).then(U_stage{ 0});
		_y >>= flow::cue_f(0x1C).then(T_event{ 1});

	//	_y >>= flow::cue_f(0x18).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x28).then(U_stage{-1});
		_y >>= flow::cue_f(0x28).then(T_event{ 1});
	//	_y >>= flow::cue_f(0x38).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x38).then(T_event{ 0});

		echo_("\nvectrol: patch");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(0 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08, 0x18);

	//	_y >>= flow::cue_f(0x08).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x08).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(U_stage{ 1});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 1});

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(1 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08);

	}
	/***/
}


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
