#pragma once
#include "./any.cc"
#include "./gate.hh"
#include "./staged.hh"
#include "./prewarped.hh"


#include "./filter.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("filter")
{
	TRY_("instantiation")
	{
		using _fit = bond::fit<>;
		using meta = any<filter<>>;

		using U_alpha = typename _fit::alpha_type;
		//\
		using SVF = confined_t<staged<0>, filter<>>;
		using SVF = prewarped_t<ordinal_constant_t<1>
		,	typename meta::refade_type::template attend<>
		,	staged< 0>
		,	filter  <>
		>;

		//\
		using Z = processor::monomer_t<prewarped<ordinal_constant_t<1>>, SVF>;
		using Z = processor::monomer_t<SVF>;

		SVF svf{};
		svf <<= occur::resample_f(44100);
		svf <<= typename  meta::  limit_type{0};
		svf <<= typename  meta::  order_type{2};
		svf <<= typename  meta::  patch_type{0};
		svf <<= typename  meta::refade_type{0};
	
		U_alpha constexpr omega = 2*2*3*3*5*5*7;
		U_alpha constexpr   rho = 1;
		U_alpha constexpr    up = 1;
		U_alpha constexpr    dn = 0;

		U_alpha _LP0{};
		U_alpha _LP1{};

		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;

		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;

		_std::array<U_alpha, 4> u_{1, 2, 0, 1};
		_std::array<U_alpha, 4> f_{3, 3, 3, 3};

		//\
		auto z = Z::bind_f(u_, f_);
//		auto z = Z::bind_f(processor::let_f(u_), processor::let_f(f_));
//		TRUE_(true);

	}
}
/***/

////////////////////////////////////////////////////////////////////////////////

TAG_("filter-ring")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;

	using U_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	U_alpha constexpr omega = 2*2*2*3*5*5*7;
	U_alpha constexpr   rho = 1;
	U_alpha constexpr    up = 1;
	U_alpha constexpr    dn = 0;

	/**/
	TRY_("filter-ring monophony")
	{
		using A_filter = filter<>;
		using A = any<A_filter>;

		using   limit_type = typename A::   limit_type;
		using   order_type = typename A::   order_type;
		using   patch_type = typename A::   patch_type;
		using   stage_type = typename A::   stage_type;
		using   shape_type = typename A::   shape_type;
		using  refade_type = typename A::  refade_type;
		using  redamp_type = typename A::  redamp_type;
		using reshape_type = typename A:: reshape_type;

		//\
		using Z_packet = stage_type;
		using Z_packet = flow::packet_t<stage_type, redamp_type>;
		using Z_event  = flow::cue_s<Z_packet>;

		flow::packed_t<          > constexpr note{};
		flow::packed_t<stage_type> constexpr note_on { 0};
		flow::packed_t<stage_type> constexpr note_off{ 1};
		flow::packed_t<stage_type> constexpr note_out{-1};

		using Z_process = prewarped_t<ordinal_constant_t<0>, gate<-1>
		,	typename A::redamp_type::template attend<>
		,	typename A::refade_type::template attend<>
	//	,	ring    <>
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

		z <<=  limit_type{0};
		z <<=  order_type{2};
		z <<=  patch_type{0};
		z <<= redamp_type{1};
		z <<= refade_type{1};

		z <<= z_sample;
		z <<= z_resize;

		z >>= occur::stage_f(-1);
		/*/
		z >>= flow::cue_f(0x08).then((note_on , redamp_type{0}));
		z >>= flow::cue_f(0x10).then((note_on , redamp_type{0}));
		z >>= flow::cue_f(0x27).then((note_off, redamp_type{root_f<-2>(2.F)}));
		/*/
		z >>= flow::cue_f(0x08).then(Z_packet{ 0, 0.0F});
		z >>= flow::cue_f(0x10).then(Z_packet{ 0, 0.0F});
		//\
		z >>= flow::cue_f(0x20).then(Z_packet{ 1, root_f<-2>(2.F)});
		z >>= flow::cue_f(0x20).then(Z_packet{ 1, 1.0F});
		/***/

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(0 == z.influx(occur::stage_f( 1)));

		echo_plot_<28>(z.store(), 0x08, 0x10);

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(1 == z.influx(occur::stage_f( 1)));

		echo_plot_<28>(z.store(), 0x00);

	}
	/***/
	/**/
	TRY_("filter-ring polyphony")
	{
		using A_filter = filter<>;
		using A = any<A_filter>;

		using   limit_type = typename A::   limit_type;
		using   order_type = typename A::   order_type;
		using   patch_type = typename A::   patch_type;
		using   stage_type = typename A::   stage_type;
		using   shape_type = typename A::   shape_type;
		using  refade_type = typename A::  refade_type;
		using  redamp_type = typename A::  redamp_type;
		using reshape_type = typename A:: reshape_type;
		using U_alpha = typename bond::fit<>::alpha_type;

		using U_payload = flow::key_s<U_stage  >;
		using U_event   = flow::cue_s<U_payload>;
		
		using U1_cue  = flow::cue_s<flow::cue_s<>>;

		using Z_process = prewarped_t<ordinal_constant_t<0>, gate<-1>
		,	typename A::redamp_type::template attend<>
		,	typename A::refade_type::template attend<>
		//\
		,	ring    <>
		,	stage_type::template assignment<redamp_type>
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

		z <<=  limit_type{0};
		z <<=  order_type{2};
		z <<=  patch_type{0};
		z <<= redamp_type{1};
		z <<= refade_type{1};

		z <<= flow::assign_f(stage_type{ 0}) << redamp_type{         (0.F)};
		z <<= flow::assign_f(stage_type{ 1}) << redamp_type{         (1.F)};
		z <<= flow::assign_f(stage_type{-1}) << redamp_type{root_f<2>(2.F)};

		z <<= z_sample;
		z <<= z_resize;

		auto const up1 = U_payload(1,  0);
		auto const dn1 = U_payload(1, -1);

		z.lead() >>= U_stage(-1);
	//	z >>= flow::cue_f(0x08) << U_payload(1,  0);
	//	z >>= flow::cue_f(0x18) << U_payload(1,  0);
		z >>= flow::cue_f(0x08).then(U_payload{1,  0});
		z >>= flow::cue_f(0x18).then(U_payload{1,  0});
	//	z >>= flow::cue_f(0x28).then(U_payload{1,  1});// Inlined below...
		z >>= flow::cue_f(0x40).then(U_payload{1,  0});
		z >>= flow::cue_f(0x48).then(U_payload{1, -1});

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		{
			echo_plot_<28>(z.store(), 0x08, 0x18);

		//	TRUE_(2 >= z.ensemble().size());// Still decaying...
		}
		z >>= flow::cue_f(0x08).then(U_payload{1,  1});
		TRUE_(0 == z.efflux(z_cursor++));
		{
			echo_plot_<28>(z.store(), 0x08);

		//	TRUE_(2 >= z.ensemble().size());// Still decaying...
		}
	//	z >>= flow::cue_f(0x00).then(U_payload{1,  0});
	//	z >>= flow::cue_f(0x08).then(U_payload{1, -1});
		TRUE_(0 == z.efflux(z_cursor++));
		{
			echo_plot_<28>(z.store(), 0x08);

		//	TRUE_(1 >= z.ensemble().size());// Still decaying...
		}
		TRUE_(0 == z.efflux(z_cursor++));
		{
		//	echo_plot_<28>(z.store());

		//	TRUE_(0 == z.ensemble().size());
		}

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
