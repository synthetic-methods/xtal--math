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
		using E_etc = process::traits_t<E_def>;
		using E_pkt = typename E_etc::stage_type;
		using E_prx = confined_t<void
		,	prewarped<_0>, gate<0>
		,	staged<-1>
		,	staged< 0>
		,	typename E_etc::redamp_type::template   attend<>
		,	typename E_etc::refade_type::template   attend<>
		,	typename E_etc::             template   attach<>
		,	typename E_etc::             template dispatch<>
		,	vectrol<>
		,	E_def
		>;
		using E_pxr = processor::monomer_t<E_prx
		,	Z_slice::template accept<E_pkt>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = E_pxr::bind_f(processor::let_f(e_omega));
		z <<= typename E_etc::  order_type{2};
		z <<= typename E_etc:: redamp_type{1};
		z <<= typename E_etc:: refade_type{0.5};
		z <<= typename E_etc::reshape_type{typename E_etc::shape_type{0.125, one - 0.125}};

		z <<= z_sample;
		z <<= z_resize;
		z >>= typename E_etc::stage_type{-1};

		z >>= flow::cue_f(0x08).then(E_pkt{ 0});
		z >>= flow::cue_f(0x18).then(E_pkt{ 0});
		z >>= flow::cue_f(0x28).then(E_pkt{-1});
	//	z >>= flow::cue_f(0x38).then(E_pkt{ 0});

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
	TRY_("vectrol: multiplex")
	{
		//\
		using E_def = filter<>;
		using E_def = filter<U_alpha[2], union ENV>;
		using E_etc = process::traits_t<E_def>;
		using E_pkt = typename E_etc::stage_type;
		using E_prx = confined_t<void
	//	,	process::lift<W_alpha>
		,	prewarped<_0>, gate<0>
		,	staged<-1>
		,	staged< 0>
		,	typename E_etc::redamp_type::template   attend<>
	//	,	typename E_etc::refade_type::template   attend<>
		,	typename E_etc::             template   attach<>
		,	typename E_etc::             template dispatch<>
		,	vectrol<>
		,	E_def
		>;
		//\
		using E_pxr = processor::conferred_t<E_prx
		using E_pxr = processor::monomer_t<E_prx
	//	,	Z_slice::template accept<E_pkt>
		,	provision::stored  <unit_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		using O_def = filter<U_alpha[2], union RING>;
		using O_etc = process::traits_t<O_def>;
		using O_pkt = flow::packet_t<typename O_etc::stage_type, typename O_etc::redamp_type>;
		//\
		using O_prx = prewarped_t<_0, gate<-1>
		using O_prx = process::confined_t<void
		,	multiplied<>
		,	prewarped<_0, gate<-1>>
	//	,	staged<-1>
		,	staged< 0>
		,	typename O_etc::redamp_type::template   attend<>
		,	typename O_etc::refade_type::template   attend<>
		,	typename O_etc::             template   attach<>
		,	typename O_etc::             template dispatch<>
		,	O_def
		>;
		//\
		using Q_vtx = atom::brace_t<U_alpha, U_alpha>;
		using Q_vtx = bond::pack_t<U_alpha, W_alpha>;
		using Q_mtx = atom::brace_t<Q_vtx[2]>;
		using Q_prx = multiplex_t<O_prx>::template matrix_t<Q_mtx>;
		using Q_pxr = processor::monomer_t<Q_prx
		,	Z_slice::template accept<O_pkt>
		,	provision::stored  <null_type[0x100]>
	//	,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr r_omega =   2*3*3*5*5*7;
		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto _1 = processor::let_f(U_alpha{one});
		auto _e = processor::let_f(e_omega);
		auto _r = processor::let_f(r_omega);
		auto _x = E_pxr::bind_f(_e);
		//\
		auto _y = Q_pxr::bind_f(_1, _x);
		auto _y = Q_pxr::bind_f(_1, E_pxr::bind_f(_e));

		//\
		_y <<= occur::math::indent_s<Q_mtx>({r_omega, 1.0});
	//	_y <<= occur::math::indent_s<Q_mtx>({r_omega, W_alpha{111., 111.}});
		_y <<= occur::math::indent_s<Q_mtx, 1>({r_omega, W_alpha{1111, 1111}});
		_y <<= occur::math::indent_s<Q_mtx, 0>({0.0    , W_alpha{0.707, 0.707}});

		_y <<= typename E_etc::  order_type{2};
		_y <<= typename E_etc:: redamp_type{1};
		_y <<= typename E_etc:: refade_type{0.5};
		_y <<= typename E_etc::reshape_type{typename E_etc::shape_type{0.125, one - 0.825}};

		_y <<= typename O_etc::  order_type{2};
		_y <<= typename O_etc:: redamp_type{1};
		_y <<= typename O_etc:: refade_type{1};

		_y <<= z_sample;
		_y <<= z_resize;
	//	_y >>= typename E_etc::stage_type{-1};
		_y >>= typename O_etc::stage_type{-1};

		_y >>= flow::cue_f(0x08).then(O_pkt{ 0});
	//	_y >>= flow::cue_f(0x08).then(E_pkt{ 0});
	//	_y >>= flow::cue_f(0x18).then(O_pkt{ 1});
	//	_y >>= flow::cue_f(0x18).then(E_pkt{ 1});

		_y >>= flow::cue_f(0x18).then(O_pkt{ 0});
	//	_y >>= flow::cue_f(0x18).then(E_pkt{ 0});
		_y >>= flow::cue_f(0x1C).then(O_pkt{ 1});

	//	_y >>= flow::cue_f(0x18).then(E_pkt{ 0});
	//	_y >>= flow::cue_f(0x18).then(O_pkt{ 0});
	//	_y >>= flow::cue_f(0x28).then(E_pkt{-1});
		_y >>= flow::cue_f(0x28).then(O_pkt{ 1});
	//	_y >>= flow::cue_f(0x38).then(E_pkt{ 0});
	//	_y >>= flow::cue_f(0x38).then(O_pkt{ 0});

		echo_rule_<28>("\u2500");

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(0 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08, 0x18);

	//	_y >>= flow::cue_f(0x08).then(O_pkt{ 0});
	//	_y >>= flow::cue_f(0x08).then(E_pkt{ 0});
	//	_y >>= flow::cue_f(0x18).then(E_pkt{ 1});
	//	_y >>= flow::cue_f(0x18).then(O_pkt{ 1});

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(1 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08);

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
