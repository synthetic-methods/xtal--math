	//\
	using V_shape  =  atom::couple_t<V_alpha[2]>;
	using V_shape  =  atom::quantity_plus_multiplies_t<V_alpha[2]>;
	using V_coeff  =  atom::math::               dot_t<V_alpha[2]>;

	using U_speed  =  occur::resample_t<>;
	using U_stage  =  occur::stage_t<>;
	using U_shape  =  occur::reinferred_t<V_shape, union SHAPE>;
	using U_coeff  =  occur::reinferred_t<V_coeff, union COEFF>;

	using X_note_scalar =  flow::packed_t<V_alpha, V_alpha, V_alpha, V_alpha, union NOTE>;
	using X_note_vector =  flow::packed_t<V_shape, V_shape, V_shape, V_shape, union NOTE>;
	using X_note_matrix =  atom::bracket_t<X_note_vector[1]>;



	using S_phasor = bond::compose<void
	,	process::phasor<V_alpha[2]>
	,	typename U_speed::attach<>
	>;
	using T_phasor = processor::monomer_t<void// Mapped and patched?
//	,	process::phasor_t<V_alpha[2], typename U_speed::template attach<>>
	,	process::confined_t<void
		,	U_frequency::affix<>
		,	S_phasor
		>
	>;

	using S_note_gate = process::zavalishin::filter<V_alpha[2], union ENV>;
	using U_note_gate = process::occurrence_t<S_note_gate>;
	using T_note_gate = processor::monomer_t<void
	,	process::confined_t<flow::mask<3>
		,	patch_t<X_note_matrix> ::   refit  <>
		,	patch_t<X_note_matrix> ::     fit  <U_shape>
		,	reuse<0, -1>
		,	coefficient_t<U_coeff> ::   attach <>
		,	pulse_t<1>             ::   infix  <>
		,	U_shape                ::   affix  <>
	//	,	U_stage                ::   attach <>
		,	U_note_gate            ::   attach <>
		,	U_note_gate            :: dispatch <>
		,	vactrol<1>
		,	S_note_gate
		,	scheme::math::zavalishin::distorted<identity>
		>
//	,	scheme::stored  <null_type[0x100]>
//	,	scheme::spooled <null_type[0x100]>
	>
	::	bind_t<T_phasor>;

	using S_note_trig = process::zavalishin::filter<V_alpha[2], union ENV>;
	using U_note_trig = process::occurrence_t<S_note_trig>;
	using T_note_trig = processor::monomer_t<void
	,	process::confined_t<flow::mask<5>
		,	patch_t<X_note_matrix> ::   refit  <>
		,	patch_t<X_note_matrix> ::     fit  <U_shape>
		,	reuse<0>
		,	coefficient_t<U_coeff> ::   attach <>
		,	pulse_t<0>             ::   infix  <>
		,	U_shape                ::   affix  <>
		,	U_stage                ::   attach <>
		,	U_note_trig            ::   attach <>
		,	U_note_trig            :: dispatch <>
		,	vactrol<0>
		,	S_note_trig
		,	scheme::math::zavalishin::distorted<identity>
		>
//	,	scheme::stored  <null_type[0x100]>
//	,	scheme::spooled <null_type[0x100]>
	>
	::	bind_t<T_phasor>;


	using V_pitch = V_alpha;
	using V_roll  = V_alpha;
	using V_lift  = V_alpha;
	using V_skew  = V_alpha;
	using V_tilt  = V_alpha;

	using X_wave_vector  =  flow:: packed_t<V_pitch, V_roll, V_lift, V_skew, V_tilt>;
	using X_wave_matrix  =  atom::bracket_t<X_wave_vector[2]>;

	using T_wave = processor::monomer_t<void
	:	process::confined_t<void
		,	codeficient<U_cross>   ::  attach  < >
		,	coefficient<U_coeff>   ::  attach  < >
		,	patch_t<X_wave_matrix> ::   refix  < >
		,	left<S_phasor>
		,	process::meth          ::    zinc  <2>
		,	process::math::zavalishin::filter_t<U_aphex[2]>   ::   defix  < >
		,	process::meth          ::  zigzag  <1>
		>
	,	typename Z_occurs :: dispatch<>// Bind at the  block-level...
	>
	::	template bind_t<T_note_gate, T_note_trig>;

	using Y_wave = processor::polymer_t<T_wave
	//\
	,	scheme::spooled<extent_constant_t<  -1>>
	,	scheme::spooled<extent_constant_t<0x20>>
	,	Z_slice::template suspend<X_note_scalar>
	>;


		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method_body(
			atom::math::phason_simplex_q auto phi_, simplex_field_q auto pitch
		,	atom::math::phason_simplex_q auto psi_, simplex_field_q auto lift
		,	auto &&...oo
		)

		template <int N_dup=0, int N_ect=0, int N_sel=0, int N_der=0>
		requires in_v<N_dup, 1> and in_v<N_der, 1>
		XTAL_DEF_(return,inline,let)
		method_body(atom::math::phason_simplex_q auto t_
		,	complex_field_q auto &&t
		,	atom::math::pade::uniplex_q auto &&s_
		,	simplex_variable_q auto f_skew
		,	simplex_variable_q auto f_tilt
		)

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(atom::math::phason_q auto &&t_, auto &&x
		method(auto &&x, atom::math::phason_simplex_q auto &&t_
		,	complex_field_q auto const &s
		,	auto &&...oo
		)	const


	z_wave <<= flow::cue_f(11).then(X_note_scalar{101, 0.5, 0.5, 0.5, 0.5});
	z_wave <<= flow::cue_f(11) << flow::mark_f(101) << X_note_scalar{0.5, 0.5, 0.5, 0.5};






	using T_wave = processor::polymer_t<void
	,	scheme::spooled <null_type[0x100]>
	>;
	confined_t<Px_wave> px_zeta; px_zeta <<= U_speed{44100};



	U_aphex const t{0, 0};
	U_aphex const s{0, 3};

	X_patched::bind_t<T_note_gate, T_note_trig>;
