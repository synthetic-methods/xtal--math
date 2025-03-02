#pragma once
#include "./any.cc"
#include "./gate.hh"
#include "./staged.hh"
#include "./prewarped.hh"


#include "./vectrol.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("vectrol")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using W_alpha = atom::math::dot_t<U_alpha[2]>;
	using Z_slice = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("vectrol: monophony")
	{
		//\
		using E_def = filter<U_alpha[2], union ENV>;
		using E_def = filter<>;
		using E_env = any_t<E_def>;
		using E_eve = typename E_env::stage_type;
		using E_pro = confined_t<void
		,	prewarped<_0>, gate<0>
		,	staged<-1>
		,	staged< 0>
		,	typename E_env::redamp_type::template   attend<>
		,	typename E_env::refade_type::template   attend<>
		,	typename E_env::             template   attach<>
		,	typename E_env::             template dispatch<>
		,	vectrol<>
		,	E_def
		>;
		using E_prx = processor::monomer_t<E_pro
		,	Z_slice::template inqueue<E_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = E_prx::bind_f(processor::let_f(e_omega));
		z <<= typename E_env::  limit_type{0};
		z <<= typename E_env::  order_type{2};
		z <<= typename E_env::  patch_type{0};
		z <<= typename E_env:: redamp_type{1};
		z <<= typename E_env:: refade_type{0.5};
		z <<= typename E_env::reshape_type{typename E_env::shape_type{0.125, one - 0.125}};

		z <<= z_sample;
		z <<= z_resize;
		z >>= typename E_env::stage_type{-1};

		z >>= flow::cue_f(0x08).then(E_eve{ 0});
		z >>= flow::cue_f(0x18).then(E_eve{ 0});
		z >>= flow::cue_f(0x28).then(E_eve{-1});
	//	z >>= flow::cue_f(0x38).then(E_eve{ 0});

		echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08, 0x18);

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(1 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08);

	}
	/***/
	/**/
	TRY_("vectrol: cross")
	{
		//\
		using E_def = filter<>;
		using E_def = filter<U_alpha[2], union ENV>;
		using E_env = any_t<E_def>;
		using E_eve = typename E_env::stage_type;
		using E_pro = confined_t<void
	//	,	process::lift<W_alpha>
		,	prewarped<_0>, gate<0>
		,	staged<-1>
		,	staged< 0>
		,	typename E_env::redamp_type::template   attend<>
	//	,	typename E_env::refade_type::template   attend<>
		,	typename E_env::             template   attach<>
		,	typename E_env::             template dispatch<>
		,	vectrol<>
		,	E_def
		>;
		using E_prx = processor::monomer_t<E_pro
		,	Z_slice::template inqueue<E_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

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
		,	Z_slice::template inqueue<R_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 2*2*2*3*5*5*7;
		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		//\
		using Q_vtx = atom::quanta_t<U_alpha, U_alpha>;
		using Q_vtx = atom::quanta_t<U_alpha, W_alpha>;
		using Q_mtx = atom::quanta_t<Q_vtx[1]>;
		using Q_pro = cross_t<Q_mtx, R_pro>;
	//	using Q_prx = cross_t<Q_mtx, R_prx>;
		using Q_prx = processor::monomer_t<Q_pro
		,	Z_slice::template inqueue<R_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		auto _1 = processor::let_f(U_alpha{one});
		auto _e = processor::let_f(e_omega);
		auto _r = processor::let_f(r_omega);
		auto _x = E_prx::bind_f(_e);
		auto _y = Q_prx::bind_f(_1, _x);

		//\
		_y <<= occur::math::indent_s<Q_mtx>({r_omega, 1});
		_y <<= occur::math::indent_s<Q_mtx>({r_omega, W_alpha{0.5, 0.5}});

		_y <<= typename E_env::  limit_type{0};
		_y <<= typename E_env::  order_type{2};
		_y <<= typename E_env::  patch_type{0};
		_y <<= typename E_env:: redamp_type{1};
		_y <<= typename E_env:: refade_type{0.5};
		_y <<= typename E_env::reshape_type{typename E_env::shape_type{0.125, one - 0.125}};

		_y <<= typename R_env::  limit_type{0};
		_y <<= typename R_env::  order_type{2};
		_y <<= typename R_env::  patch_type{0};
		_y <<= typename R_env:: redamp_type{1};
		_y <<= typename R_env:: refade_type{1};

		_y <<= z_sample;
		_y <<= z_resize;
		_y >>= typename E_env::stage_type{-1};
		_y >>= typename R_env::stage_type{-1};

		_y >>= flow::cue_f(0x08).then(E_eve{ 0});
		_y >>= flow::cue_f(0x08).then(R_eve{ 0});
		_y >>= flow::cue_f(0x18).then(E_eve{ 0});
		_y >>= flow::cue_f(0x18).then(R_eve{ 0});
		_y >>= flow::cue_f(0x28).then(E_eve{-1});
		_y >>= flow::cue_f(0x28).then(R_eve{-1});
	//	_y >>= flow::cue_f(0x38).then(E_eve{ 0});

		echo_rule_<28>("\u2500");

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(0 == _y.influx(occur::stage_f(-1)));

		echo_plot_<28>(_y.store(), 0x08, 0x18);
//
//		TRUE_(0 == _y.efflux(z_cursor++));
//	//	TRUE_(1 == _y.influx(occur::stage_f(-1)));
//
//		echo_plot_<28>(_y.store(), 0x08);

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
