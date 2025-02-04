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
	TRY_("instantiation")
	{
		using _fix = bond::fixture<>;
		using U_alpha = typename _fix::alpha_type;

		using U_stage = occur::stage_t<>;
		using U_chunk = schedule::chunk_t<provision::spooled<extent_constant_t<0x10>>>;
		using U_value = flow::packet_t<U_stage>;
		using U_event = flow::cue_s<U_value>;

		//\
		using SVF = prewarped_t<ordinal_constant_t<0>, ring<>>;
		using SVF = confined_t<prewarped<ordinal_constant_t<0>>, ring<>>;

		using Z = processor::monomer_t<SVF, provision::stored<>>;

		SVF svf{};
		svf <<= typename occur::sample_t<>{44100};
		svf <<= typename ring<>::    limit_type{0};
		svf <<= typename ring<>::    order_type{2};
		svf <<= typename ring<>::   select_type{0};
		svf <<= typename ring<>:: topology_type{0};
		svf <<= typename ring<>::  damping_type{1};
		svf <<= typename ring<>::  balance_type{1};
	
		U_alpha constexpr omega = 2*2*2*3*5*5*7;
		U_alpha constexpr   rho = 1;
		U_alpha constexpr    up = 1;
		U_alpha constexpr    dn = 0;

		U_alpha _LP0{};
		U_alpha _LP1{};

		_std::array<U_alpha, 0x100> u_; u_.fill(up);// u_[0] = up;
		_std::array<U_alpha, 0x100> f_; f_.fill(omega);

		//\
		auto z = Z::bind_f(u_, f_);
		auto z = Z::bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x20);
		auto z_render = occur::render_t<>(0x20);
//
		z <<= typename occur::sample_t<>{44100};
		z <<= typename ring<>::    limit_type{0};
		z <<= typename ring<>::    order_type{2};
		z <<= typename ring<>::   select_type{0};
		z <<= typename ring<>:: topology_type{0};
		z <<= typename ring<>::  damping_type{1};
		z <<= typename ring<>::  balance_type{1};

		z <<= z_resize;

		z <<= occur::stage_f(0);

		z >>= z_render++;

		plot<5>(z.store());

	}
	TRY_("instantiation")
	{
		using _fix = bond::fixture<>;
		using U_alpha = typename _fix::alpha_type;

		using U_chunk = schedule::chunk_t<provision::spooled<extent_constant_t<0x10>>>;

		using U_stage = occur::stage_t<>;
		using U_value = flow::key_s<U_stage>;
		using U_event = flow::cue_s<U_value>;

		//\
		using Z_process = prewarped_t<ordinal_constant_t<0>, ring<>>;
		using Z_process = confined_t<prewarped<ordinal_constant_t<0>>, ring<>>;

		using Z_processor = processor::polymer_t<Z_process
		,	U_chunk::template inqueue<U_value>
		,	provision::stored<>
		,	provision::spooled<>
		>;

		U_alpha constexpr omega = 2*2*3*3*5*5*7;

		_std::array<U_alpha, 0x100> f_; f_.fill(omega);
//		Z_processor::template bind_t<> z_{processor::let_f(f_)};
//
//		z_ <<= U_resize(8);
//
//		z_ <<= U_stage(-1);
//
//		z_ <<= U_event(2, 65, 1);
	
	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
