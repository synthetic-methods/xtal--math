#pragma once
#include "./any.hh"

#include "../../provision/saturated.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Integrating filter based on the Topology Preserving Transform (TPT), \
described in _The Art of VA Synthesis_ by Vadim Zavalishin. \

///\
The non-linearity is supplied as a `process`-`template` using `provision::saturated`. \
The process must conform to the signature `<M_ism, M_car>`, \
defining a stateless `method_f<N_var, ...>` within the `subtype`. \

///\
The parameters `M_ism` and `M_car` determine the type and return-value of curve, respectively. \
`M_ism` is expected to yield convex/concave curves for positive/negative values, \
and `M_car` is expected to return the slope when `== -1`. \

///\note\
Despite the parameterization defined by `any<filter<...>>`, \
`filter<...>::method` is polymorphic and can accomodate anything up to the given cache-width \
(determined by `U_pole[N_pole][2]`). \
\
For example, with base-types of `double` and `std::complex<double>` respectively, \
the storage required is `16` and `32` bytes-per-pole. \

template <class ..._s>	struct  filter;
template <class ..._s>	concept filter_q = bond::tag_p<filter, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class U_pole, auto N_pole>
struct any<filter<U_pole[N_pole]>>
{
//	USED
	using        scope_type =  atom::couple_t<typename bond::fit<U_pole>::alpha_type[N_pole]>;
	using        state_type =  atom::couple_t<U_pole[N_pole]>;
	using        stage_type = occur::stage_t<>;

//	ATTACHED
	using        input_type = occur::inferred_t<struct   INPUT, valued_u<state_type>>;
	using        scale_type = occur::inferred_t<struct   SCALE, valued_u<scope_type>>;
	using      damping_type = occur::inferred_t<struct DAMPING, valued_u<scope_type>>;
	using      balance_type = occur::inferred_t<struct BALANCE, valued_u<scope_type>>;
	using      rescope_type = occur::inferred_t<struct RESCOPE,          scope_type >;
	using         zoom_type = occur::inferred_t<struct    ZOOM, valued_u<scope_type>>;
	using     sampling_type = occur::sampling_t<>;

//	DISPATCHED
	using       select_type = occur::inferred_t<struct  SELECT, unsigned int, bond::word<  2>>;
	//\
	using        order_type = occur::inferred_t<struct   ORDER, unsigned int, bond::seek_s<1 + state_type::size()>>;
	using        order_type = occur::inferred_t<struct   ORDER, unsigned int, bond::word  <1 + state_type::size()>>;
	using        patch_type = occur::inferred_t<struct   PATCH, unsigned int, bond::word  <2>>;
	using        limit_type = occur::inferred_t<struct   LIMIT, unsigned int, bond::word  <2>>;

	using superkind = provision::voiced<void
	,	typename   select_type::template dispatch<>
	,	typename    order_type::template dispatch<>
	,	typename patch_type::template dispatch<>
	,	typename    limit_type::template dispatch<>// TODO: Replace with `shape`'s parameters?
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	//	USED
		using        scope_type = any::       scope_type;
		using        state_type = any::       state_type;
		using        stage_type = any::       stage_type;

	//	ATTACHED
		using        input_type = any::       input_type;
		using        scale_type = any::       scale_type;
		using      damping_type = any::     damping_type;
		using      balance_type = any::     balance_type;
		using         zoom_type = any::        zoom_type;
		using     sampling_type = any::    sampling_type;

	//	DISPATCHED
		using       select_type = any::      select_type;
		using        order_type = any::       order_type;
		using     patch_type = any::    patch_type;
		using        limit_type = any::       limit_type;


	};
};
template <scalar_q A>
struct any<filter<A>>
:	any<filter<A[2]>>
{
};
template <>
struct any<filter<>>
:	any<filter<typename bond::fit<>::alpha_type>>
{
};


////////////////////////////////////////////////////////////////////////////////

template <vector_q A>
struct filter<A>
{
	using          metakind = any<filter<A>>;
	using        state_type = typename metakind::        state_type;
	using        stage_type = typename metakind::        stage_type;
	using     sampling_type = typename metakind::     sampling_type;
	using         zoom_type = typename metakind::         zoom_type;

	using superkind = bond::compose<bond::tag<filter>
	,	provision::memorized<state_type, state_type>
	,	metakind
	,	typename    stage_type::template expect<>
	,	typename sampling_type::template attach<>
	,	typename     zoom_type::template attach<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <int N_sel=0, int N_ord=0, int N_pat=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x_input
		,	absolve_u<decltype(x_input)> s_scale
		,	absolve_u<decltype(x_input)> s_damping
		,	absolve_u<decltype(x_input)> y_balance=0
		)
		noexcept -> auto
		{
			using X = XTAL_ALL_(x_input);
			using W = absolve_u<X>;

			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return XTAL_REF_(x_input);
			}
			XTAL_0IF (0 == N_pat) {
				using U_outputs = atom::couple_t<X[N_ord + 1]>;
				using U_scalars = atom::couple_t<W[N_ord + 1]>;
				union U_io {U_outputs outputs; U_scalars scalars;};

				U_io io{[=] () XTAL_0FN -> U_scalars {
					auto const  &u = cut_f<[] XTAL_1FN_(value) (bond::fit<X>::minilon_f(0))>(s_damping);
					auto const u02 =             two*u , u04 = u02*two;
					auto const u12 = term_f(one, two,u), w24 = u02*u12;
					XTAL_IF0
					XTAL_0IF (1 == N_ord) {return {one, one};}
					XTAL_0IF (2 == N_ord) {return {one, u02, one};}
					XTAL_0IF (3 == N_ord) {return {one, u12, u12, one};}
					XTAL_0IF (4 == N_ord) {return {one, u04, w24, u04, one};}
				}()};

				(void) edit_f<N_ord, N_pat, Ns...>(XTAL_REF_(x_input), s_scale, io);// io.scalars -> io.outputs

				auto constexpr I_ =  static_cast<unsigned>(N_sel);
				auto constexpr I0 = _std::countr_one(I_ >>  0) +  0, J0 = I0 + 1;
				auto constexpr I1 = _std::countr_one(I_ >> J0) + J0, J1 = I1 + 1;
				W const &y1 = y_balance;
				//\
				W const  y0 = one;
				W const  y0 = term_f<-1, 2>(one, y1);

				if constexpr (I_ == 0) {
					static_assert(I0 == 0);
					static_assert(I1 == 1);
				}
				XTAL_IF0
				XTAL_0IF (0  == N_ord) {return x_input;}
				XTAL_0IF (I1 >= N_ord) {return io.outputs.self(constant_t<1>{}).product(atom::couple_t<W[1]>{y0    });}
				XTAL_0IF (I1 <  N_ord) {return io.outputs.self(constant_t<2>{}).product(atom::couple_t<W[2]>{y0, y1});}
			}
			XTAL_0IF (1 == N_pat) {
				return XTAL_REF_(x_input);
			}
		}
		template <int N_ord=0, int N_pat=0, int N_lim=0, auto ...Ns> requires (1 <= N_ord and N_pat == 0)
		XTAL_DEF_(inline,let)
		edit_f(auto const &x_input
		,	absolve_u<decltype(x_input)> s_scale
		,	auto &io
		)
		noexcept -> void
		{
			using X = XTAL_ALL_(x_input);
			using W = absolve_u<X>;

			using U_scalars_ = atom::couple_t<W[N_ord]>;
			using U_outputs_ = atom::couple_t<X[N_ord]>;
			using U_poles_   = atom::couple_t<X[N_ord]>;

			auto  scalars = io.scalars; auto scalars_ = scalars.self(cardinal_constant_t<N_ord>{});
			auto &outputs = io.outputs; auto outputs_ = outputs.self(cardinal_constant_t<N_ord>{});

			auto  memorized  = S_::template memory<U_poles_, U_poles_>();
		//	(void) cut_t<[] XTAL_1FN_(value) (-bond::fit<X>::diplo_f(7))>::edit_f(memorized);
			auto &states_ = get<0>(memorized);
			auto &slopes_ = get<1>(memorized);

		//	Initialize `slopes*`:
			slopes_.template blanket<1>();
			slopes_ *= scalars_;

		//	Initialize `outputs*`:
			get<0 >(outputs) = get<0>(slopes_);
			bond::seek_forward_f<N_ord - 1>([&] (auto I) XTAL_0FN {
				get<I + 1>(outputs_) = term_f(get<I + 1>(slopes_), get<I>(outputs_), s_scale);
			});
			get<N_ord>(outputs) = root_f<-1, (1)>(term_f(one, get<N_ord - 1>(outputs_), s_scale));

			auto constexpr K_lim = provision::saturated_q<S_> and (0 < N_lim);

		//	Integrate `states*` with `scalars*` and `outputs*`:
			bond::seek_forward_f<1 + 1*K_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {get<N_ord>(outputs) *= x_input - outputs_.product(states_);}
				XTAL_0IF (1 <= K) {get<N_ord>(outputs)  = x_input - outputs_.product(slopes_);}

				bond::seek_backward_f<N_ord>([&] (auto I) XTAL_0FN {
					XTAL_IF0
					XTAL_0IF (0 == K_lim) {
						get<I>(outputs_) = term_f(get<I>(states_), get<I + 1>(outputs), s_scale);
					}
					XTAL_0IF (1 == K_lim) {
						auto const exput = term_f(get<I>(states_), get<I + 1>(outputs), s_scale);
						auto const slope = S_::template saturate_t<-1, -1>::template method_f<N_lim, Ns...>(exput) + (get<I>(scalars_) - one);
						get<I>(outputs_) = exput;
						get<I>( slopes_) = slope;
					}
				});
			});

		//	Finalize states and `outputs*`/`scalars*`:
			bond::seek_forward_f<N_ord>([&] (auto I) XTAL_0FN {
				get<I>(states_) = term_f(-get<I>(states_), two, get<I>(outputs_));
			});
			slopes_ /= scalars_;
		//	outputs_ *= slopes_;// TODO: Make this line optional?
		}

	};
};
template <scalar_q A>
struct filter<A>
:	filter<A[2]>
{
};
template <>
struct filter<>
:	filter<typename bond::fit<>::alpha_type>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
