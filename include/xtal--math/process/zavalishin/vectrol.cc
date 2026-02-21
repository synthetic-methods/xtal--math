#pragma once
#include "./any.cc"
#include "./pulse.hh"
#include "./reuse.hh"
#include "../../occur/all.hh"
#include "../../provision/prewarping.hh"

#include "./vectrol.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test::XTAL_NUM
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
using namespace _test;

struct reshape
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		using shape_type      = typename S_::shape_type;
	//	using shape_parameter = typename S_::shape_parameter;
		using shape_parameter = occur::math::produce_s<atom::math::dot_t<shape_type[2][2]>, union SHAPE>;
		using shape_parametry = occur::math::indent_s<shape_parameter, 0>;
		using shape_occurence = occur::math::indent_s<shape_parameter, 1>;

	};
};


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
		using S_content = filter<>;
		using S_content = filter<U_alpha[2], union ENV>;
		//\
		using S_context = confined_t<S_content>;
		using S_context = occur::context_t<S_content>;
		using S_stage   = typename S_context::stage_type;
		using S_shape   = typename S_context::shape_type;
		using S_process = confined_t<void
		,	provision::math::prewarping< 0>
		,	pulse< 0>
		,	reuse< 0>
		,	reuse<-1>
		,	occur::math::zavalishin::dash_t<S_shape>::template attend<>
	//	,	vectrol<>// TODO: Annotate as below...
	//	,	typename S_context::damp_parameter::template   attend<>
		,	typename S_context::fade_parameter::template   attend<>
		,	typename S_context::                template   attach<>
		,	typename S_context::                template dispatch<>
	//	,	reshape
		,	S_content
		>;
		using S_processor = processor::monomer_t<S_process
		,	Z_slice::template accept<S_stage>
		,	provision::stored  <null_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		U_alpha constexpr e_omega = 2*2*3*3*5*5;
		auto z_resize = occur::resize_t<>(0x020);
		auto z_cursor = occur::cursor_t<>(0x020);
		auto z_sample = occur::resample_f(44100);

		auto z = S_processor::bind_f(processor::let_f(e_omega));
		z <<= typename S_context::order_attribute{2};
	//	z <<= typename S_context:: damp_parameter{1};
		z <<= typename S_context:: fade_parameter{0.5};
	//	z <<= typename S_context::shape_parameter{-0. - one, 0.};
		z <<= occur::math::zavalishin::dash_t<S_shape>{one - 0.0, zero - 1.0};

		z <<= z_sample;
		z <<= z_resize;
		z >>= typename S_context::stage_type{-1};

		z >>= flow::cue_f(0x08).then(S_stage{ 0});
		z >>= flow::cue_f(0x18).then(S_stage{ 0});
		z >>= flow::cue_f(0x28).then(S_stage{-1});
	//	z >>= flow::cue_f(0x38).then(S_stage{ 0});

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
	/*/
	TRY_("vectrol: patch")
	{
		//\
		using S_content = filter<>;
		using S_content = filter<U_alpha[2], union ENV>;
		using S_context = occur::context_t<S_content>;
		using S_stage   = typename S_context::stage_type;
		using S_process = confined_t<void
	//	,	process::lift<W_alpha>
		,	provision::math::prewarping< 0>
		,	pulse< 0>
		,	reuse< 0>
		,	reuse<-1>
		,	vectrol<>
	//	,	typename S_context::damp_parameter::template   attend<>
	//	,	typename S_context::fade_parameter::template   attend<>
		,	typename S_context::                template   attach<>
		,	typename S_context::                template dispatch<>
		,	S_content
		>;
		//\
		using S_processor = processor::conferred_t<S_process
		using S_processor = processor::monomer_t<S_process
	//	,	Z_slice::template accept<S_stage>
		,	provision::stored  <unit_type[0x100]>
		,	provision::spooled <null_type[0x100]>
		>;

		using T_content = filter<U_alpha[2], union RING>;
		using T_context = occur::context_t<T_content>;
		using T_stage   = typename T_context::stage_type;
		using T_event   = flow::packet_t<T_stage, typename T_context::damp_parameter>;

		using T_process = process::confined_t<void
		,	coefficient<>
		,	provision::math::prewarping< 0>
		,	pulse<-1>
		,	reuse< 0>
	//	,	reuse<-1>
		,	typename T_context::damp_parameter::template   attend<>
		,	typename T_context::fade_parameter::template   attend<>
		,	typename T_context::                template   attach<>
		,	typename T_context::                template dispatch<>
		,	T_content
		>;
		//\
		using X_vector    = atom::brace_t<U_alpha, U_alpha>;
		using X_vector    = bond::pack_t<U_alpha, W_alpha>;
		using X_matrix    = atom::brace_t<X_vector[2]>;
		using X_process   = patch_t<T_process>::template matrix_t<X_matrix>;
		using X_processor = processor::monomer_t<X_process
		,	Z_slice::template accept<T_event>
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
		_y <<= occur::math::indent_s<X_matrix, 0>({0.0    , W_alpha{0.707, 0.707}});

		_y <<= typename S_context::order_attribute{2};
//		_y <<= typename S_context:: damp_parameter{1};
		_y <<= typename S_context:: fade_parameter{0.5};
		_y <<= typename S_context::shape_parameter{0.0 - one,  1.0};

		_y <<= typename T_context::order_attribute{2};
		_y <<= typename T_context:: damp_parameter{1};
		_y <<= typename T_context:: fade_parameter{1};

		_y <<= z_sample;
		_y <<= z_resize;
	//	_y >>= typename S_context::stage_type{-1};
		_y >>= typename T_context::stage_type{-1};

		_y >>= flow::cue_f(0x08).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x08).then(S_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 1});
	//	_y >>= flow::cue_f(0x18).then(S_stage{ 1});

		_y >>= flow::cue_f(0x18).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x18).then(S_stage{ 0});
		_y >>= flow::cue_f(0x1C).then(T_event{ 1});

	//	_y >>= flow::cue_f(0x18).then(S_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x28).then(S_stage{-1});
		_y >>= flow::cue_f(0x28).then(T_event{ 1});
	//	_y >>= flow::cue_f(0x38).then(S_stage{ 0});
	//	_y >>= flow::cue_f(0x38).then(T_event{ 0});

		echo_("\nvectrol: patch");
	//	echo_rule_<28>("\u2500");

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(0 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08, 0x18);

	//	_y >>= flow::cue_f(0x08).then(T_event{ 0});
	//	_y >>= flow::cue_f(0x08).then(S_stage{ 0});
	//	_y >>= flow::cue_f(0x18).then(S_stage{ 1});
	//	_y >>= flow::cue_f(0x18).then(T_event{ 1});

		TRUE_(0 == _y.efflux(z_cursor++));
	//	TRUE_(1 == _y.influx(occur::stage_f(-1)));
		echo_plot_<28>(_y.store(), 0x08);

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
