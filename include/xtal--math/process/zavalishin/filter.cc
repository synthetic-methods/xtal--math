#pragma once
#include "./any.cc"
#include "../../occur/all.hh"
#include "../../provision/all.hh"
#include "./reuse.hh"


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

	using U_resample = occur::resample_t<>;
	using U_resync   = occur::resync_t<>;


	using U_coeff  = atom::math::dot_t<U_alpha[2]>;
	using X_coeff  = occur::inferred_t<U_coeff, union COEFF>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("instantiation")
	{
		using R_def = filter<>;
		using R_etc = occur::codex_t<R_def>;
		using R_prx = confined_t<void
		,	per_t<U_resample>          ::   refix <1>
		,	U_resync                   ::   attach <>
		,	coefficient_t<X_coeff>     ::   attach <>
	//	,	reuse< 0>
		,	R_etc                      ::   attach <>
		,	R_etc                      :: dispatch <>
		,	R_def
		,	provision::math::zavalishin::shaped<identity>
		>;
		using R_pxr = processor::monomer_t<R_prx>;

		R_prx svf{};
		svf <<= occur::resample_f(44100);
		svf <<= typename R_etc::  order_attribute{2};
		svf <<= X_coeff{1, 0};
	
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
		,	typename per_t<U_resample>::template refix<1>
		//\
		,	filter<union RAMP>
		,	filter<U_alpha[2], union RAMP>
		>;
		using Y_ring = confined_t<void
		,	typename per_t<U_resample>::template refix<1>
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

	using U_resample = occur::resample_t<>;
	using U_resync = occur::resync_t<>;
	using U_stage    = occur::stage_t<>;

	using U_coeff  = atom::math::dot_t<U_alpha[2]>;
	using X_coeff  = occur::inferred_t<U_coeff, union COEFF>;

	using Y_trig  = pulse_t< 0>;
	using Y_gate  = pulse_t< 1>;
	using Y_hold  = pulse_t<-1>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("filter-ring monophony")
	{
		using R_def = filter<U_alpha[2], union RING>;
		using R_etc = occur::codex_t<R_def>;
		using R_eve = flow::packet_t<U_stage, typename R_etc::damp_parameter>;
		using R_prx = confined_t<void
		,	per_t<U_resample>::refix<0>
		//\
		,	reuse<0, -1>
		,	reuse<   -1>
		,	U_resync               ::   attach <>
		,	coefficient_t<X_coeff> ::   attach <>
		,	Y_hold                 ::   infix  <>
		,	R_etc::damp_parameter  ::   affix  <>
		,	R_etc::                     attach <>
		,	R_etc::                   dispatch <>
		,	R_def
		>;
		using R_pxr = processor::monomer_t<R_prx
		,	Z_slice::template suspend<R_eve>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 2*2*2*3*5*5*7;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = R_pxr::bind_f(processor::let_f(r_omega));
		z <<= typename R_etc::order_attribute{2};
		z <<= typename R_etc:: damp_parameter{1};
		z <<= X_coeff{0, 1};

		z <<= typename occur::resync_t<>{0};

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
		TRUE_(1 == z.influx(occur::stage_f( 0)));// Would be unchanged...
		TRUE_(0 == z.influx(occur::stage_f( 1)));
		TRUE_(0 == z.influx(occur::stage_f(-1)));
		echo_plot_<28>(z.store(), 0x08, 0x10);

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.influx(occur::stage_f( 0)));
		TRUE_(1 == z.influx(occur::stage_f( 1)));// Would be unchanged...
		TRUE_(0 == z.influx(occur::stage_f(-1)));
		//\
		echo_plot_<28>(z.store(), 0x00);
		echo_plot_<28>(z.store(), 0x01);

	}
	/***/
	/**/
	TRY_("filter-ring polyphony")
	{
		using R_def = filter<U_alpha[2], union RING>;
		using R_etc = occur::codex_t<R_def>;
		using R_eve = flow::key_s<U_stage>;

		using R_prx = confined_t<void
	//	,	reuse<0, -1>
		,	reuse<   -1>
		,	per_t<U_resample>      ::   refix <0>
		,	coefficient_t<X_coeff> ::   attach <>
		,	U_resync               ::   attach <>
		,	Y_gate                 ::   infix  <>
		,	U_stage::template assignment<typename R_etc::damp_parameter>
		,	R_etc::damp_parameter  ::   affix  <>
		,	R_etc                  ::   attach <>
		,	R_etc                  :: dispatch <>
		,	R_def
		>;
		using R_pxr = processor::polymer_t<R_prx
		,	Z_slice::template suspend<R_eve>
		,	provision::stored <null_type[0x100]>
		,	provision::spooled<null_type[0x100]>
		>;

		U_alpha constexpr r_omega = 59*61;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = R_pxr::bind_f(processor::let_f(r_omega));
		z <<= typename R_etc::order_attribute{2};
		z <<= typename R_etc:: damp_parameter{1};
		z <<= X_coeff{0, 1};
		z <<= typename occur::resync_t<>{0};
		z <<= flow::assign_f(U_stage{ 0}) << typename R_etc::damp_parameter{0.000F};
		z <<= flow::assign_f(U_stage{ 1}) << typename R_etc::damp_parameter{0.060F};
		z <<= flow::assign_f(U_stage{-1}) << typename R_etc::damp_parameter{0.707F};

		z <<= z_sample;
		z <<= z_resize;

		z.lead() >>= U_stage{-1};

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


////////////////////////////////////////////////////////////////////////////////

TAG_("vectrol")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	using W_alpha = atom::math::dot_t<U_alpha[2]>;
	using Z_slice = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

	using U_resample = occur::resample_t<>;
	using U_stage = occur::stage_t<>;

	using U_coeff = atom::math::dot_t<U_alpha[2]>;
	using X_coeff = occur::inferred_t<U_coeff, union COEFF>;

	using Y_trig  = pulse_t< 0>;
	using Y_gate  = pulse_t< 1>;
	using Y_hold  = pulse_t<-1>;

	using _0 = ordinal_constant_t<0>;
	using _1 = ordinal_constant_t<1>;

	/**/
	TRY_("vectrol: monophony")
	{
		//\
		using S_content = filter<>;
		using S_content = filter<U_alpha[2], union ENV>;
		//\
		using S_codex = confined_t<S_content>;
		using S_codex = occur::codex_t<S_content>;

		using S_damp_   = occur::math::zavalishin::probe_t<typename S_codex::codata_type>;
	//	using S_damp    = typename S_codex::damp_parameter;

		using S_order   = typename S_codex::order_attribute;

		using S_process = confined_t<void
		,	reuse<0, -1>
		,	coefficient_t<X_coeff> ::   attach <>
		,	per_t<U_resample>      ::   refix <0>
		,	Y_trig                 ::   infix  <>
		,	S_damp_                ::   affix  <>
	//	,	S_damp                 ::   affix  <>
		,	S_codex                ::   attach <>
		,	S_codex                :: dispatch <>
		,	S_content
		>;
		using S_processor = processor::monomer_t<S_process
		,	Z_slice::template suspend<U_stage>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = S_processor::bind_f(processor::let_f(e_omega));
		z <<= S_order{2};
		z <<= S_damp_{one - 0.0, zero - 1.0};
		z <<= X_coeff{1, 1};

		z <<= z_sample;
		z <<= z_resize;
		z >>= U_stage{-1};

		z >>= flow::cue_f(0x08).then(U_stage{ 0});
		z >>= flow::cue_f(0x18).then(U_stage{ 0});
		z >>= flow::cue_f(0x28).then(U_stage{-1});
	//	z >>= flow::cue_f(0x38).then(U_stage{ 0});

		echo_("\nvectrol: monophony");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == z.efflux(z_cursor++));
		TRUE_(0 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08, 0x18);

		TRUE_(0 == z.efflux(z_cursor++));
	//	TRUE_(1 == z.influx(occur::stage_f(-1)));

		echo_plot_<28>(z.store(), 0x08);

	}
	/***/
	/**/
	TRY_("vectrol: patch")
	{
		using S_content = filter<U_alpha[2], union ENV>;
		using S_codex   = occur::codex_t<S_content>;
		using S_damp_   = occur::math::zavalishin::probe_t<typename S_codex::codata_type>;
		using S_damp    = typename S_codex::  damp_parameter;
		using S_order   = typename S_codex:: order_attribute;

		using S_process = confined_t<void
		,	reuse< 0, -1>
	//	,	process::lift<W_alpha>
		,	typename per_t<U_resample> ::   refix <0>
		,	typename Y_trig            ::   infix  <>
		,	typename S_damp_           ::   affix  <>
		,	typename S_codex           ::   attach <>
		,	typename S_codex           :: dispatch <>
		,	S_content
		>;
		//\
		using S_processor = processor::conferred_t<S_process
		using S_processor = processor::monomer_t<S_process
	//	,	Z_slice::template suspend<U_stage>
		,	provision::stored  <unit_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		using T_content = filter<U_alpha[2], union RING>;
		using T_codex   = occur::codex_t<T_content>;
	//	using T_damp_   = occur::math::zavalishin::probe_t<typename T_codex::codata_type>;
		using T_damp    = typename T_codex:: damp_parameter;

		using T_order   = typename T_codex::order_attribute;
		using T_event   = flow::packet_t<U_stage, T_damp>;

		using T_process = process::confined_t<void
		,	reuse< 0>
		,	coefficient_t<X_coeff> ::   attach <>
		,	coefficient_t<       > ::   unfix  <>
		,	per_t<U_resample>      ::   refix <0>
		,	Y_gate                 ::   infix  <>
	//	,	T_damp_                ::   affix  <>
		,	T_damp                 ::   affix  <>
		,	T_codex                ::   attach <>
		,	T_codex                :: dispatch <>
		,	T_content
		>;
		using X_vector    = atom::brace_t<U_alpha, W_alpha>;
		using X_matrix    = atom::brace_t<X_vector[2]>;
		using X_process   = patch_t<T_process>::template matrix_t<X_matrix>;
		using X_processor = processor::monomer_t<X_process
		,	Z_slice::template suspend<T_event>
		,	provision::stored  <null_type[0x100]>
	//	,	provision::spooled <null_type[0x100]>
		>;

		static_assert(           fungible_q<typename occur::math::indent_s<X_matrix, 1>::data_type, X_matrix>);
		static_assert(occur::math::indent_q<typename occur::math::indent_s<X_matrix, 1>           , X_matrix>);

		U_alpha constexpr r_omega = 3*3*3*5*5*5;
		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto _1 = processor::let_f(U_alpha{one});
		auto _e = processor::let_f(e_omega);
		auto _r = processor::let_f(r_omega);
		auto _x = S_processor::bind_f(_e);
		//\
		auto _y = X_processor::bind_f(_1, _x);
		auto _y = X_processor::bind_f(_1, S_processor::bind_f(_e));

		//\
		_y <<= occur::math::indent_s<X_matrix>({r_omega, 1.0});
	//	_y <<= occur::math::indent_s<X_matrix>({r_omega, W_alpha{111., 111.}});
		_y <<= occur::math::indent_s<X_matrix, 1>({r_omega, W_alpha{1111, 1111}});
		_y <<= occur::math::indent_s<X_matrix, 0>({0.0    , W_alpha{1.00, 1.00}});

		_y <<= S_order{2};
		_y <<= S_damp_{one - 0.0, zero - 1.0};
	//	_y <<= S_damp {1.0};

		_y <<= T_order{2};
	//	_y <<= T_damp_{1};
		_y <<= T_damp {1};
		_y <<= X_coeff {0, 1};

		_y <<= z_sample;
		_y <<= z_resize;
		_y >>= U_stage{-1};

		_y >>= flow::cue_f(0x08).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x08).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 1});
	//	_y >>= flow::cue_f(0x18).then(U_stage{ 1});

		_y >>= flow::cue_f(0x18).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x18).then(U_stage{ 0});
		_y >>= flow::cue_f(0x1C).then(T_event{ 1});

	//	_y >>= flow::cue_f(0x18).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x28).then(U_stage{-1});
		_y >>= flow::cue_f(0x28).then(T_event{ 1});
	//	_y >>= flow::cue_f(0x38).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x38).then(T_event{ 0});

		echo_("\nvectrol: patch");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(0 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08, 0x18);

	//	_y >>= flow::cue_f(0x08).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x08).then(U_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(U_stage{ 1});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 1});

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(1 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08);

	}
	/***/
}


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
