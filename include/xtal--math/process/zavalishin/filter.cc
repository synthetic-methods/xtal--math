#pragma once
#include "./any.cc"
#include "./intake.hh"
#include "./retake.hh"
#include "../../provision/prewarping.hh"


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
TAG_("filter")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using W_alpha = atom::math::dot_t<U_alpha[2]>;
	using Z_slice = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("instantiation")
	{
		using R_def = filter<>;
		using R_etc = occur::context_t<R_def>;
		using R_prx = confined_t<void
		,	provision::math::prewarping< 1>
		,	retake< 0>
		,	typename R_etc::fade_parameter::template   attend<>
	//	,	typename R_etc::zoom_parameter::template   attend<>
		,	typename R_etc::             template   attach<>
		,	typename R_etc::             template dispatch<>
		,	R_def
		,	provision::math::saturation<identity>
		>;
		using R_pxr = processor::monomer_t<R_prx>;

		R_prx svf{};
		svf <<= occur::resample_f(44100);
		svf <<= typename R_etc::  order_attribute{2};
		svf <<= typename R_etc:: fade_parameter{0};
	
		U_alpha constexpr r_omega = 2*2*3*3*5*5*7;
		U_alpha constexpr   rho = 1;
		U_alpha constexpr    up = 1;
		U_alpha constexpr    dn = 0;

		U_alpha _LP0{};
		U_alpha _LP1{};

		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, r_omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;

		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, r_omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;

		_std::array<U_alpha, 4> u_{1, 2, 0, 1};
		_std::array<U_alpha, 4> f_{3, 3, 3, 3};

	}
	TRY_("filter parameterization")
	{
		using Y_ramp = confined_t<void
		,	provision::math::prewarping< 1>
		//\
		,	filter<union RAMP>
		,	filter<U_alpha[2], union RAMP>
		>;
		using Y_ring = confined_t<void
		,	provision::math::prewarping< 1>
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
	using Z_slice = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("filter-ring monophony")
	{
		using R_def = filter<U_alpha[2], union RING>;
		using R_etc = occur::context_t<R_def>;
		using R_eve = flow::packet_t<typename R_etc::stage_type, typename R_etc::damp_parameter>;
		using R_prx = confined_t<void
		,	provision::math::prewarping< 0>
		,	intake<-1>
		,	retake< 0>
		,	retake<-1>
		,	typename R_etc::damp_parameter::template   attend<>
		,	typename R_etc::fade_parameter::template   attend<>
		,	typename R_etc::             template   attach<>
		,	typename R_etc::             template dispatch<>
		,	R_def
		>;
		using R_pxr = processor::monomer_t<R_prx
		,	Z_slice::template accept<R_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 2*2*2*3*5*5*7;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = R_pxr::bind_f(processor::let_f(r_omega));
		z <<= typename R_etc::  order_attribute{2};
		z <<= typename R_etc:: damp_parameter{1};
		z <<= typename R_etc:: fade_parameter{1};

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
		using R_def = filter<U_alpha[2], union RING>;
		using R_etc = occur::context_t<R_def>;
		using R_eve = flow::key_s<typename R_etc::stage_type>;
		using R_prx = confined_t<void
		,	provision::math::prewarping< 0>
		,	intake< 1>
		,	retake< 0>
		,	retake<-1>
		,	typename R_etc:: stage_type::template assignment<typename R_etc::damp_parameter>
		,	typename R_etc::damp_parameter::template   attend<>
		,	typename R_etc::fade_parameter::template   attend<>
		,	typename R_etc             ::template   attach<>
		,	typename R_etc             ::template dispatch<>
		,	R_def
		>;
		using R_pxr = processor::polymer_t<R_prx
		,	Z_slice::template accept<R_eve>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 59*61;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = R_pxr::bind_f(processor::let_f(r_omega));
		z <<= typename R_etc::  order_attribute{2};
		z <<= typename R_etc:: damp_parameter{1};
		z <<= typename R_etc:: fade_parameter{1};
		z <<= flow::assign_f(typename R_etc::stage_type{ 0}) << typename R_etc::damp_parameter{0.000F};
		z <<= flow::assign_f(typename R_etc::stage_type{ 1}) << typename R_etc::damp_parameter{0.060F};
		z <<= flow::assign_f(typename R_etc::stage_type{-1}) << typename R_etc::damp_parameter{0.707F};

		z <<= z_sample;
		z <<= z_resize;

		z.lead() >>= typename R_etc::stage_type{-1};

		z >>= flow::cue_f(0x08).then(R_eve{69, 0});
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
		echo_plot_<28>(z.store(), 0x08, 0x18);

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


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
