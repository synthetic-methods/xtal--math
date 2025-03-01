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
The parameters `M_ism` and `M_car` determine the type and return-value of shape, respectively. \
`M_ism` is expected to yield convex/concave shapes for positive/negative values, \
and `M_car` is expected to return the slope when `== -1`. \

///\note\
Despite the parameterization defined by `any<filter<...>>`, \
`filter<...>::method` is polymorphic and can accomodate anything up to the given cache-width \
(determined by `U_pole[N_pole][2]`). \
\
For example, with base-types of `double` and `std::complex<double>` respectively, \
the storage required is `16` and `32` bytes-per-pole. \

template <class ..._s>	struct  filter;
template <class ..._s>	concept filter_q = bond::tagged_with_p<filter, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct any<filter<_s...>>
{
	using superkind = any<class_template<_s...>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
	
	public:
		using S_::S_;

		template <extent_type N_mask=-1>
		struct   attach
		{
			template <class R>
			using subtype = bond::compose_s<R
			,	typename T_::    stage_type::template  inspect<N_mask>
			,	typename T_::   rezoom_type::template   attach<N_mask>
			,	typename T_:: resample_type::template   attach<N_mask>
			>;

		};
		template <extent_type N_mask=-1>
		struct dispatch
		{
			template <class R>
			using subtype = bond::compose_s<R, provision::voiced<void
			,	typename T_::    order_type::template dispatch<N_mask>
			,	typename T_::    patch_type::template dispatch<N_mask>
			,	typename T_::    limit_type::template dispatch<N_mask>// TODO: Replace with `shape`'s parameters?
			>>;

		};

	};
};
template <scalar_q A>
struct any<filter<A>> : any<filter<A[2]>>
{
};
template <>
struct any<filter< >> : any<filter<typename bond::fit<>::alpha_type>>
{
};


////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct filter
{
	using metakind = any  <filter>;
	using metatype = any_t<filter>;

	using state_type = typename metatype::state_type;
	using slope_type = typename metatype::slope_type;

	using superkind = bond::compose<bond::tag<filter>
	,	metakind
	,	provision::memorized<state_type, slope_type>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <int N_ord=0, int N_pat=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_u<decltype(x)> s_gain
		,	unstruct_u<decltype(x)> s_damp
		,	unstruct_u<decltype(x)> y_fade
		)
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x); using X2 = atom::couple_t<X[2]>;
			using W = unstruct_u<X>; using W2 = atom::couple_t<W[2]>;

			return dot_f(method<N_ord, N_pat, Ns...>(XTAL_REF_(x), s_gain, s_damp)
			,	W2{term_f<-1, 2>(one, y_fade), y_fade}
			);//TODO: Replace with `fade`.
		}
		template <int N_ord=0, int N_pat=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_u<decltype(x)> s_gain
		,	unstruct_u<decltype(x)> s_damp
		)
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x); using X2 = atom::couple_t<X[2]>;
			using W = unstruct_u<X>; using W2 = atom::couple_t<W[2]>;

			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return X2{XTAL_REF_(x)};
			}
			XTAL_0IF (0 == N_pat) {
				return method<N_ord, N_pat, Ns...>(XTAL_REF_(x), s_gain
				,	[=] () XTAL_0FN -> atom::couple_t<W[N_ord + 1]> {
						auto const  &u = cut_f<[] XTAL_1FN_(to) (bond::fit<X>::minilon_f(0))>(s_damp);
						auto const u02 =             two*u , u04 = u02*two;
						auto const u12 = term_f(one, two,u), w24 = u02*u12;
						XTAL_IF0
						XTAL_0IF (1 == N_ord) {return {one, one};}
						XTAL_0IF (2 == N_ord) {return {one, u02, one};}
						XTAL_0IF (3 == N_ord) {return {one, u12, u12, one};}
						XTAL_0IF (4 == N_ord) {return {one, u04, w24, u04, one};}
					}()
				);
			}
			XTAL_0IF (1 == N_pat) {
				return X2{XTAL_REF_(x)};
			}
		}
		template <int N_ord=0, int N_pat=0, int N_lim=0, auto ...Ns> requires (1 <= N_ord and N_pat == 0)
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_u<decltype(x)> s_gain
		,	atom::couple_q<unit_type[N_ord + 1]> auto &&scalars
		)
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x); using X2 = atom::couple_t<X[2]>;
			using W = unstruct_u<X>; using W2 = atom::couple_t<W[2]>;

			using U_state_ = atom::couple_t<X[N_ord]>;
			using U_slope_ = atom::couple_t<X[N_ord]>;

			atom::couple_t<X[N_ord + 1]> outputs{};

			auto outputs_ = outputs.self(cardinal_constant_t<N_ord>{});
			auto scalars_ = scalars.self(cardinal_constant_t<N_ord>{});

			auto     mem  = S_::template memory<U_state_, U_slope_>();
		//	(void) cut_t<[] XTAL_1FN_(to) (-bond::fit<X>::diplo_f(7))>::edit_f(mem);
			auto &states_ = get<0>(mem);
			auto &slopes_ = get<1>(mem);

		//	Initialize `slopes*`:
			slopes_.template blanket<1>();
			slopes_ *= scalars_;

		//	Initialize `outputs*`:
			get<0 >(outputs) = get<0>(slopes_);
			bond::seek_out_f<N_ord - 1>([&] (auto I) XTAL_0FN {
				get<I + 1>(outputs_) = term_f(get<I + 1>(slopes_), get<I>(outputs_), s_gain);
			});
			get<N_ord>(outputs) = root_f<-1, (1)>(term_f(one, get<N_ord - 1>(outputs_), s_gain));

			auto constexpr K_lim = provision::saturated_q<S_> and (0 < N_lim);

		//	Integrate `states*` with `scalars*` and `outputs*`:
			bond::seek_out_f<1 + 1*K_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {get<N_ord>(outputs) *= x - dot_f(outputs_, states_);}
				XTAL_0IF (1 <= K) {get<N_ord>(outputs)  = x - dot_f(outputs_, slopes_);}

				bond::seek_out_f<-N_ord>([&] (auto I) XTAL_0FN {
					XTAL_IF0
					XTAL_0IF (0 == K_lim) {
						get<I>(outputs_) = term_f(get<I>(states_), get<I + 1>(outputs), s_gain);
					}
					XTAL_0IF (1 == K_lim) {
						auto const exput = term_f(get<I>(states_), get<I + 1>(outputs), s_gain);
						auto const slope = S_::template saturate_t<-1, -1>::
							template method_f<N_lim, Ns...>(exput) + (get<I>(scalars_) - one);
						get<I>(outputs_) = exput;
						get<I>( slopes_) = slope;
					}
				});
			});

		//	Finalize states and `outputs*`/`scalars*`:
			bond::seek_out_f<N_ord>([&] (auto I) XTAL_0FN {
				get<I>(states_) = term_f(-get<I>(states_), two, get<I>(outputs_));
			});
			slopes_ /= scalars_;
		//	outputs_ *= slopes_;// TODO: Make this line optional?

			//\
			return X2{get<0>(outputs), get<1>(outputs)};
			return X2{outputs.self(constant_t<2>{})};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
