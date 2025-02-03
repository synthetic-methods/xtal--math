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
		using U = typename _fix::alpha_type;

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
	
		U constexpr omega = 2*2*2*3*5*5*7;
		U constexpr   rho = 1;
		U constexpr    up = 1;
		U constexpr    dn = 0;

		U _LP0{};
		U _LP1{};

		_std::array<U, 0x100> u_; u_.fill(up);// u_[0] = up;
		_std::array<U, 0x100> f_; f_.fill(omega);

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
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
