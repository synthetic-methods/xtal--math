#pragma once
#include "./any.cc"

#include "./prewarped.hh"
#include "./staged.hh"
#include "./iota.hh"

#include "./ring.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("ring")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;
	using U0_cue  = flow::cue_s<>;

	using U_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	U_alpha constexpr omega = 2*2*2*3*5*5*7;
	U_alpha constexpr   rho = 1;
	U_alpha constexpr    up = 1;
	U_alpha constexpr    dn = 0;

	/**/
	TRY_("ring: monophony")
	{
		using Z_filter = filter<>;
		using Z = any<Z_filter>;

		using refade_type = typename Z::refade_type;
		using redamp_type = typename Z::redamp_type;
		using reshape_type = typename Z::reshape_type;
		using   shape_type = typename Z::  shape_type;
		using   stage_type = typename Z::  stage_type;

		//\
		using Z_packet = stage_type;
		using Z_packet = flow::packet_t<stage_type, redamp_type>;
		using Z_event  = flow::cue_s<Z_packet>;
		using Z_cue    = flow::cue_s<>;

		using A_filter = filter<>;

		using A = any<A_filter>;
		using Z_process = prewarped_t<ordinal_constant_t<0>, iota<-1>
		,	typename A::redamp_type::template attend<>
		,	typename A::refade_type::template attend<>
		,	ring    <>
		,	staged<-1>
		,	staged< 0>
		,	A_filter
		>;
		using Z_processor = processor::monomer_t<Z_process
		,	U_slicer::template inqueue<Z_packet>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= typename A::   limit_type{0};
		z <<= typename A::   order_type{2};
		z <<= typename A::   patch_type{0};
		z <<= typename A:: redamp_type{1};
		z <<= typename A:: refade_type{1};

		z <<= z_sample;
		z <<= z_resize;
		z <<= U_stage(-1);

		z <<= Z_cue(0x08).then(Z_packet{ 0, 0});
		z <<= Z_cue(0x10).then(Z_packet{ 0, 0});
		z <<= Z_cue(0x27).then(Z_packet{-1, root_f<-2>(2.F)});

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.efflux(occur::stage_f(-1)));

		echo_rule_<25>('=');
		echo_plot_<25>(z.store());

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(1 == z.efflux(occur::stage_f(-1)));

		echo_rule_<25>('-');
		echo_plot_<25>(z.store());

	}
	/***/
	/**/
	TRY_("ring: polyphony")
	{
		using U_alpha = typename bond::fit<>::alpha_type;

		using U_payload = flow::key_s<U_stage  >;
		using U_event   = flow::cue_s<U_payload>;
		
		using U0_cue  = flow::cue_s<>;
		using U1_cue  = flow::cue_s<flow::cue_s<>>;

		using A_filter = filter<>;

		using A = any<A_filter>;
		using Z_process = prewarped_t<ordinal_constant_t<0>, iota<-1>
		,	typename A::redamp_type::template attend<>
		,	typename A::refade_type::template attend<>
		,	ring    <>
		,	staged<-1>
		,	staged< 0>
		,	filter  <>
		>;
		using Z_processor = processor::polymer_t<Z_process
		,	U_slicer::template inqueue<U_payload>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::template bind_f(processor::let_f(f_));

		auto z_resize = occur::    resize_t<>(0x020);
		auto z_cursor = occur::    cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= typename A::   limit_type{0};
		z <<= typename A::   order_type{2};
		z <<= typename A::   patch_type{0};
		z <<= typename A:: redamp_type{1};
		z <<= typename A:: refade_type{1};

		z <<= z_sample;
		z <<= z_resize;
		//\
		z <<= U_stage(-1);
		z.lead() <<= U_stage(-1);

		auto const up1 = U_payload(1,  0);
		auto const dn1 = U_payload(1, -1);

	//	z <<= U0_cue(0x08) << U_payload(1,  0);
	//	z <<= U0_cue(0x18) << U_payload(1,  0);
		z <<= U0_cue(0x08).then(U_payload{1,  0});
		z <<= U0_cue(0x18).then(U_payload{1,  0});
	//	z <<= U0_cue(0x28).then(U_payload{1, -1});
		z <<= U0_cue(0x40).then(U_payload{1,  0});
		z <<= U0_cue(0x48).then(U_payload{1, -1});


		TRUE_(0 == z.efflux(z_cursor++));
		{
			echo_rule_<25>('=');
			echo_plot_<25>(z.store());

		//	TRUE_(2 >= z.ensemble().size());// Still decaying...
		}
		z <<= U0_cue(0x08).then(U_payload{1, -1});
		TRUE_(0 == z.efflux(z_cursor++));
		{
			echo_rule_<25>('-');
			echo_plot_<25>(z.store());

		//	TRUE_(2 >= z.ensemble().size());// Still decaying...
		}
	//	z <<= U0_cue(0x00).then(U_payload{1,  0});
	//	z <<= U0_cue(0x08).then(U_payload{1, -1});
		TRUE_(0 == z.efflux(z_cursor++));
		{
			echo_rule_<25>('-');
			echo_plot_<25>(z.store());

		//	TRUE_(1 >= z.ensemble().size());// Still decaying...
		}
		TRUE_(0 == z.efflux(z_cursor++));
		{
		//	echo_rule_<25>('-');
		//	echo_plot_<25>(z.store());

		//	TRUE_(0 == z.ensemble().size());
		}

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
