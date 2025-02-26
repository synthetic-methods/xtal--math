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
/**/
TAG_("filter")
{
	TRY_("parameterization")
	{
		using T_fit = bond::fit<>;
		using T_alpha = typename T_fit::alpha_type;


		using Y_ramp = prewarped_t<ordinal_constant_t<1>
		,	filter<T_alpha[2], union RAMP>
		>;
		using Y_ring = prewarped_t<ordinal_constant_t<1>
		,	filter<T_alpha[2], union RING>
		>;
		static_assert(different_q<typename Y_ramp::order_type, typename Y_ring::order_type>);

	};
	TRY_("instantiation")
	{
		using T_fit = bond::fit<>;
		using T_alpha = typename T_fit::alpha_type;

		using A_filter = filter<T_alpha[2], union MAIN>;

		using meta = any<A_filter>;

		//\
		using SVF = confined_t<staged<0>, filter<>>;
		using SVF = prewarped_t<ordinal_constant_t<1>
		,	typename meta::refade_type::template attend<>
		,	staged< 0>
		,	A_filter
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
	
		T_alpha constexpr omega = 2*2*3*3*5*5*7;
		T_alpha constexpr   rho = 1;
		T_alpha constexpr    up = 1;
		T_alpha constexpr    dn = 0;

		T_alpha _LP0{};
		T_alpha _LP1{};

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

		_std::array<T_alpha, 4> u_{1, 2, 0, 1};
		_std::array<T_alpha, 4> f_{3, 3, 3, 3};

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
	using T_fit   = bond::fit<>;
	using T_alpha = typename T_fit::alpha_type;
	using T_delta = typename T_fit::delta_type;
	using T_sigma = typename T_fit::sigma_type;

	using U_stage = occur::stage_t<>;
	using U_key   = flow::key_s<>;

	using U_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	T_alpha constexpr omega = 2*2*2*3*5*5*7;
	T_alpha constexpr   rho = 1;
	T_alpha constexpr    up = 1;
	T_alpha constexpr    dn = 0;

	/**/
	TRY_("filter-ring monophony")
	{
		using A_filter = filter<>;
		using A = any<A_filter>;

		using   input_type = typename A::   input_type;
		using   limit_type = typename A::   limit_type;
		using   order_type = typename A::   order_type;
		using   patch_type = typename A::   patch_type;
		using   stage_type = typename A::   stage_type;
		using   curve_type = typename A::   curve_type;
		using  refade_type = typename A::  refade_type;
		using  redamp_type = typename A::  redamp_type;
		using recurve_type = typename A:: recurve_type;

		using Z_packet = flow::packet_t<stage_type, redamp_type>;
		using Z_event  = flow::cue_s<Z_packet>;

		using Z_process = prewarped_t<ordinal_constant_t<0>, gate<-1>
		,	redamp_type::template attend<>
		,	refade_type::template attend<>
		,	staged<-1>
		,	staged< 0>
		,	A_filter
		>;
		using Z_processor = processor::monomer_t<Z_process
		,	U_slicer::template inqueue<Z_packet>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<T_alpha, 0x100> f_; f_.fill(omega);
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
		z >>= flow::cue_f(0x08).then(Z_packet{ 0, 0.0F});
		z >>= flow::cue_f(0x10).then(Z_packet{ 0, 0.0F});
		//\
		z >>= flow::cue_f(0x20).then(Z_packet{ 1, root_f<-2>(2.F)});
		z >>= flow::cue_f(0x20).then(Z_packet{ 1, 1.0F});

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
		using   curve_type = typename A::   curve_type;
		using  refade_type = typename A::  refade_type;
		using  redamp_type = typename A::  redamp_type;
		using recurve_type = typename A:: recurve_type;
		using T_alpha = typename bond::fit<>::alpha_type;

		using V_event = flow::key_s<U_stage>;
		using U_event = flow::cue_s<V_event>;
		
		using U0_cue  = flow::cue_s<      >;
		using U1_cue  = flow::cue_s<U0_cue>;

		using Z_process = prewarped_t<ordinal_constant_t<0>, gate<-1>
		,	A::redamp_type::template attend<>
		,	A::refade_type::template attend<>
		,	stage_type::template assignment<redamp_type>
		,	staged<-1>
		,	staged< 0>
		,	filter  <>
		>;
		using Z_processor = processor::polymer_t<Z_process
		,	U_slicer::template inqueue<V_event>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		_std::array<T_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::template bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
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

		auto const up1 = V_event(1,  0);
		auto const dn1 = V_event(1, -1);

		z.lead() >>= U_stage(-1);

		z >>= flow::cue_f(0x08).then(V_event{1,  0});
		z >>= flow::cue_f(0x18).then(V_event{1,  0});
	//	z >>= flow::cue_f(0x28).then(V_event{1,  1});// Inlined below...
		z >>= flow::cue_f(0x40).then(V_event{1,  0});
		z >>= flow::cue_f(0x48).then(V_event{1, -1});

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.ensemble().size());
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08, 0x18);

		z >>= flow::cue_f(0x08).then(V_event{1,  1});

		TRUE_(2 >= z.ensemble().size());// Still decaying...
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08);

	//	z >>= flow::cue_f(0x00).then(V_event{1,  0});
	//	z >>= flow::cue_f(0x08).then(V_event{1, -1});

		TRUE_(2 >= z.ensemble().size());// Still decaying...
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08);

		TRUE_(1 >= z.ensemble().size());// Still decaying...
		TRUE_(0 == z.efflux(z_cursor++));
	//	echo_plot_<28>(z.store());

	//	TRUE_(0 == z.ensemble().size());

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
