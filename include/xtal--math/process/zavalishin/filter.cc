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


template <int M_ism=0, int M_car=0>
struct identishape;

template <int M_ism>
struct identishape<M_ism, -0> : bond::compose<discarded<1>, identishape<M_ism, -1>>
{
};
template <int M_ism>
struct identishape<M_ism, -1>
{
	template <class S>
	class subtype : public S
	{
	public:
		using S::S;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&x, auto &&...oo)
		noexcept -> auto
		{
			return XTAL_ALL_(x) {one};
		}

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
		using R_env = any_t<R_def>;
		using R_pro = prewarped_t<_1
		,	staged< 0>
		,	typename R_env::refade_type::template   attend<>
	//	,	typename R_env::rezoom_type::template   attend<>
		,	typename R_env::             template   attach<>
		,	typename R_env::             template dispatch<>
		,	R_def
		,	provision::saturated<identishape>
		>;
		//\
		using R_prx = processor::monomer_t<prewarped<_1>, R_pro>;
		using R_prx = processor::monomer_t<R_pro>;

		R_pro svf{};
		svf <<= occur::resample_f(44100);
		svf <<= typename R_env::  limit_type{0};
		svf <<= typename R_env::  order_type{2};
		svf <<= typename R_env::  patch_type{0};
		svf <<= typename R_env:: refade_type{0};
	
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
		using Y_ramp = prewarped_t<_1
		//\
		,	filter<union RAMP>
		,	filter<U_alpha[2], union RAMP>
		>;
		using Y_ring = prewarped_t<_1
		//\
		,	filter<union RING>
		,	filter<U_alpha[2], union RING>
		>;
		static_assert(different_q<typename Y_ramp::order_type, typename Y_ring::order_type>);

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
		using R_env = any_t<R_def>;
		using R_eve = flow::packet_t<typename R_env::stage_type, typename R_env::redamp_type>;
		using R_pro = prewarped_t<_0, gate<-1>
		,	staged<-1>
		,	staged< 0>
		,	typename R_env::redamp_type::template   attend<>
		,	typename R_env::refade_type::template   attend<>
		,	typename R_env::             template   attach<>
		,	typename R_env::             template dispatch<>
		,	R_def
		>;
		using R_prx = processor::monomer_t<R_pro
		,	Z_slice::template accept<R_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 2*2*2*3*5*5*7;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = R_prx::bind_f(processor::let_f(r_omega));
		z <<= typename R_env::  limit_type{0};
		z <<= typename R_env::  order_type{2};
		z <<= typename R_env::  patch_type{0};
		z <<= typename R_env:: redamp_type{1};
		z <<= typename R_env:: refade_type{1};

		z <<= z_sample;
		z <<= z_resize;

		z >>= occur::stage_f(-1);
		z >>= flow::cue_f(0x08).then(R_eve{ 0, 0.0F});
		z >>= flow::cue_f(0x10).then(R_eve{ 0, 0.0F});
		//\
		z >>= flow::cue_f(0x20).then(R_eve{ 1, 0.5F});
		z >>= flow::cue_f(0x21).then(R_eve{ 1, 0.5F});

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
		using R_def = filter<U_alpha[2], union RING>;
		using R_env = any_t<R_def>;
		using R_eve = flow::key_s<typename R_env::stage_type>;
		using R_pro = prewarped_t<_0, gate<1>
		,	staged< 0>
		,	staged<-1>
		,	typename R_env:: stage_type::template assignment<typename R_env::redamp_type>
		,	typename R_env::redamp_type::template   attend<>
		,	typename R_env::refade_type::template   attend<>
		,	typename R_env             ::template   attach<>
		,	typename R_env             ::template dispatch<>
		,	R_def
		>;
		using R_prx = processor::polymer_t<R_pro
		,	Z_slice::template accept<R_eve>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 2*2*2*3*5*5*7;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = R_prx::bind_f(processor::let_f(r_omega));
		z <<= typename R_env::  limit_type{0};
		z <<= typename R_env::  order_type{2};
		z <<= typename R_env::  patch_type{0};
		z <<= typename R_env:: redamp_type{1};
		z <<= typename R_env:: refade_type{1};
		z <<= flow::assign_f(typename R_env::stage_type{ 0}) << typename R_env::redamp_type{          (0.F)};
		z <<= flow::assign_f(typename R_env::stage_type{ 1}) << typename R_env::redamp_type{root_f<-1>(2.F)};
		z <<= flow::assign_f(typename R_env::stage_type{-1}) << typename R_env::redamp_type{root_f<-2>(2.F)};

		z <<= z_sample;
		z <<= z_resize;

		z.lead() >>= typename R_env::stage_type{-1};

		z >>= flow::cue_f(0x08).then(R_eve{1,  0});
		z >>= flow::cue_f(0x18).then(R_eve{1,  0});
	//	z >>= flow::cue_f(0x28).then(R_eve{1,  1});// Inlined below...
	//	z >>= flow::cue_f(0x38).then(R_eve{1, -1});
		z >>= flow::cue_f(0x40).then(R_eve{1,  0});
		z >>= flow::cue_f(0x50).then(R_eve{1,  1});

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.ensemble().size());
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08, 0x18);

		z >>= flow::cue_f(0x08).then(R_eve{1,  1});// Inlined from above...

		TRUE_(2 >= z.ensemble().size());// Still decaying...
		TRUE_(0 == z.efflux(z_cursor++));
		echo_plot_<28>(z.store(), 0x08);

	//	z >>= flow::cue_f(0x00).then(R_eve{1,  0});// Outlined above...
	//	z >>= flow::cue_f(0x08).then(R_eve{1, -1});// Outlined above...

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
