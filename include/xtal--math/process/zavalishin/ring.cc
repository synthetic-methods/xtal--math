#pragma once
#include "./any.cc"
#include "./ring.hh"// testing...

#include "./prewarped.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("ring")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;
	using U_cue   = flow::cue_s<>;

	using U_chunk = schedule::chunk_t<provision::spooled<extent_constant_t<0x10>>>;

	U_alpha constexpr omega = 2*2*2*3*5*5*7;
	U_alpha constexpr   rho = 1;
	U_alpha constexpr    up = 1;
	U_alpha constexpr    dn = 0;

	/**/
	TRY_("ring: monophony")
	{
		using U_value =             U_stage ;
		using U_event = flow::cue_s<U_value>;

		using Y = prewarped_t<ordinal_constant_t<0>, ring<>>;
		using Z = processor::monomer_t<Y
		,	U_chunk::template inqueue<U_value>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z::bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_render = occur::render_t<>(0x020);
		auto z_sample = occur::sample_t<>(44100);

		z <<= typename ring<>::    limit_type{0};
		z <<= typename ring<>::    order_type{2};
		z <<= typename ring<>::   select_type{0};
		z <<= typename ring<>:: topology_type{0};
		z <<= typename ring<>::  damping_type{1};
		z <<= typename ring<>::  balance_type{1};

		z <<= z_sample;
		z <<= z_resize;
		z <<= U_stage(     -1);

		z <<= U_event(0x04, 0);
		z <<= U_event(0x14,-1);

		z >>= z_render++;

		plot<25>(z.store());

	}
	/***/
	/*/
	TRY_("ring: polyphony")
	{
		using U_value = flow::key_s<U_stage>;
		using U_event = flow::cue_s<U_value>;

		//\
		using Z_process = prewarped_t<ordinal_constant_t<0>, ring<>>;
		using Z_process = confined_t<prewarped<ordinal_constant_t<0>>, ring<>>;


	//	using Z_processor = processor::confined_t<void
	//	,	U_chunk::template inqueue<U_key, U_stage>
	//	,	processor::polymer<Z_process
	//		,	provision::stored <null_type[0x100]>
	//		,	provision::spooled<null_type[0x100]>
	//		>
	//	>;
		using Z_processor = processor::polymer_t<void
		,	Z_process
		,	U_chunk::template inqueue<U_value>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::template bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_render = occur::render_t<>(0x020);
		auto z_sample = occur::sample_t<>(44100);

		z <<= typename ring<>::    limit_type{0};
		z <<= typename ring<>::    order_type{2};
		z <<= typename ring<>::   select_type{0};
		z <<= typename ring<>:: topology_type{0};
		z <<= typename ring<>::  damping_type{1};
		z <<= typename ring<>::  balance_type{1};

		z <<= z_sample;
		z <<= z_resize;
		z <<= U_stage(-1);

		auto const up1 = U_value(1,  0);
		auto const dn1 = U_value(1, -1);

		echo("\n<<<");
	//	z <<=                flow::key_s<U_stage>(1,  0);
		z <<= U_event(0x04, 1,  0);
		z <<= U_event(0x14, 1, -1);
	//	z <<= U_cue(0x04) << U_key(1) << U_stage( 0);
	//	z <<= U_cue(0x14) << U_key(1) << U_stage(-1);
		echo(">>>\n");

		z >>= z_render++;

		plot<25>(z.store());

	}
	/***/
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
