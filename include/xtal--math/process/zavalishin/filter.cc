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
	using T_alpha   = typename bond::fit<>::alpha_type;
	using T_aphex   = typename bond::fit<>::aphex_type;

	using A_filter  = filter<>;
	using W_filter  = any_t<A_filter>;
	using U_limit   = typename W_filter::   limit_type;
	using U_order   = typename W_filter::   order_type;
	using U_patch   = typename W_filter::   patch_type;
	using U_stage   = typename W_filter::   stage_type;
	using U_shape   = typename W_filter::   shape_type;
	using U_refade  = typename W_filter::  refade_type;
	using U_redamp  = typename W_filter::  redamp_type;
	using U_reshape = typename W_filter:: reshape_type;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	TRY_("instantiation")
	{
		//\
		using SVF = confined_t<staged<0>, filter<>>;
		using SVF = prewarped_t<_1
		,	staged< 0>
		,	typename U_refade::template   attend<>
		,	typename W_filter::template   attach<>
		,	typename W_filter::template dispatch<>
		,	A_filter
		>;

		//\
		using Z = processor::monomer_t<prewarped<_1>, SVF>;
		using Z = processor::monomer_t<SVF>;

		SVF svf{};
		svf <<= occur::resample_f(44100);
		svf <<= U_limit{0};
		svf <<= U_order{2};
		svf <<= U_patch{0};
		svf <<= U_refade{0};
	
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
	TRY_("filter parameterization")
	{
		using Y_ramp = prewarped_t<_1
		//\
		,	filter<union RAMP>
		,	filter<T_alpha[2], union RAMP>
		>;
		using Y_ring = prewarped_t<_1
		//\
		,	filter<union RING>
		,	filter<T_alpha[2], union RING>
		>;
		static_assert(different_q<typename Y_ramp::order_type, typename Y_ring::order_type>);

	};
}
/***/

////////////////////////////////////////////////////////////////////////////////

TAG_("filter-ring")
{
	using T_alpha   = typename bond::fit<>::alpha_type;
	using T_aphex   = typename bond::fit<>::aphex_type;

	using A_filter  = filter<T_alpha[2], union RING>;
	using W_filter  = any_t<A_filter>;
	using U_limit   = typename W_filter::   limit_type;
	using U_order   = typename W_filter::   order_type;
	using U_patch   = typename W_filter::   patch_type;
	using U_stage   = typename W_filter::   stage_type;
	using U_shape   = typename W_filter::   shape_type;
	using U_refade  = typename W_filter::  refade_type;
	using U_redamp  = typename W_filter::  redamp_type;
	using U_reshape = typename W_filter:: reshape_type;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	using Z_slicer = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	/**/
	TRY_("filter-ring monophony")
	{
		using U_event = flow::packet_t<U_stage, U_redamp>;
		using Z_event = flow::cue_s<U_event>;

		using Z_process = prewarped_t<_0, gate<-1>
		,	staged<-1>
		,	staged< 0>
		,	typename U_redamp::template   attend<>
		,	typename U_refade::template   attend<>
		,	typename W_filter::template   attach<>
		,	typename W_filter::template dispatch<>
		,	A_filter
		>;
		using Z_processor = processor::monomer_t<Z_process
		,	Z_slicer::template inqueue<U_event>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		T_alpha constexpr omega = 2*2*2*3*5*5*7;

		_std::array<T_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= U_limit  {0};
		z <<= U_order  {2};
		z <<= U_patch  {0};
		z <<= U_redamp {1};
		z <<= U_refade {1};

		z <<= z_sample;
		z <<= z_resize;

		z >>= occur::stage_f(-1);
		z >>= flow::cue_f(0x08).then(U_event{ 0, 0.0F});
		z >>= flow::cue_f(0x10).then(U_event{ 0, 0.0F});
		//\
		z >>= flow::cue_f(0x20).then(U_event{ 1, 0.5F});
		z >>= flow::cue_f(0x21).then(U_event{ 1, 0.5F});

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(1 == z.influx(flow::assess_f(occur::stage_f( 0))));// Would be unchanged...
		TRUE_(0 == z.influx(flow::assess_f(occur::stage_f( 1))));
		TRUE_(0 == z.influx(flow::assess_f(occur::stage_f(-1))));
		echo_plot_<28>(z.store(), 0x08, 0x10);

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.influx(flow::assess_f(occur::stage_f( 0))));
		TRUE_(1 == z.influx(flow::assess_f(occur::stage_f( 1))));// Would be unchanged...
		TRUE_(0 == z.influx(flow::assess_f(occur::stage_f(-1))));
		//\
		echo_plot_<28>(z.store(), 0x00);
		echo_plot_<28>(z.store(), 0x01);

	}
	/***/
	/**/
	TRY_("filter-ring polyphony")
	{
		using U_event = flow::key_s<U_stage>;
		using Z_event = flow::cue_s<U_event>;
		using Z_cue   = flow::cue_s<       >;
		
		using Z_process = prewarped_t<_0, gate<-1>
		,	staged<-1>
		,	staged< 0>
		,	typename U_redamp ::template   attend<>
		,	typename U_refade ::template   attend<>
		,	typename W_filter ::template   attach<>
		,	typename W_filter ::template dispatch<>
		,	typename U_stage  ::template assignment<U_redamp>
		,	A_filter
		>;
		using Z_processor = processor::polymer_t<Z_process
		,	Z_slicer::template inqueue<U_event>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		T_alpha constexpr omega = 2*2*2*3*5*5*7;

		_std::array<T_alpha, 0x100> f_; f_.fill(omega);
		auto z = Z_processor::template bind_f(processor::let_f(f_));

		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		z <<= U_limit{0};
		z <<= U_order{2};
		z <<= U_patch{0};
		z <<= U_redamp{1};
		z <<= U_refade{1};
		z <<= flow::assign_f(U_stage{ 0}) << U_redamp{          (0.F)};
		z <<= flow::assign_f(U_stage{ 1}) << U_redamp{root_f<-1>(2.F)};
		z <<= flow::assign_f(U_stage{-1}) << U_redamp{root_f<-2>(2.F)};

		z <<= z_sample;
		z <<= z_resize;

		z.lead() >>= U_stage{-1};

		z >>= Z_cue{0x08}.then(U_event{1,  0});
		z >>= Z_cue{0x18}.then(U_event{1,  0});
	//	z >>= Z_cue{0x28}.then(U_event{1,  1});// Inlined below...
		z >>= Z_cue{0x40}.then(U_event{1,  0});
		z >>= Z_cue{0x48}.then(U_event{1, -1});

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.ensemble().size());
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08, 0x18);

		z >>= Z_cue{0x08}.then(U_event{1,  1});// Inlined from above...

		TRUE_(2 >= z.ensemble().size());// Still decaying...
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08);

	//	z >>= Z_cue{0x00}.then(U_event{1,  0});// Outlined above...
	//	z >>= Z_cue{0x08}.then(U_event{1, -1});// Outlined above...

		TRUE_(2 >= z.ensemble().size());// Still decaying...
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08);

		TRUE_(1 >= z.ensemble().size());// Still decaying...
		TRUE_(0 == z.efflux(z_cursor++));
	//	echo_plot_<28>(z.store());

		TRUE_(0 == z.ensemble().size());

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
