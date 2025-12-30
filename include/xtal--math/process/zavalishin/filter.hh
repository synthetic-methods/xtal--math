#pragma once
#include "./any.hh"

#include "./scaffold.hh"
#include "../../provision/prewarping.hh"
#include "../../provision/saturation.hh"


XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Integrating filter using the Topology Preserving Transform (TPT).

Implementation outlined in _The Art of VA Synthesis_ by Vadim Zavalishin.

The non-linearity is supplied as a `process`-`template` using `provision::math::saturation`.
The process must conform to the signature `<M_ism, M_car>`,
defining a stateless `method_f<N_var, ...>` within the `subtype`.

The parameters `M_ism` and `M_car` determine the type and return-value of shape, respectively.
`M_ism` is expected to yield convex/concave shapes for positive/negative values,
and `M_car` is expected to return the slope when `== -1`.

\note    Despite the parameterization defined by `occur::context<filter<...>>`, `filter<...>::method` is
polymorphic and can accomodate scaffold up to the given cache-width (determined by `U_pole[N_pole][2]`).
For example, with base-types of `double` and `std::complex<double>` respectively,
the storage required is `16` and `32` bytes-per-pole.
*/
template <class ..._s>	struct  filter;
template <class ..._s>	concept filter_q = bond::tag_in_p<filter, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct filter
{
	using metatype = occur::context_t<filter>;

	using state_type = typename metatype::state_type;
	using slope_type = typename metatype::slope_type;

	using superkind = bond::compose<bond::tag<filter>
	,	provision::memorized<state_type, slope_type>
	,	scaffold<As...>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;

		template <int N_ism=0, int N_car=0, auto ...Ns>
		XTAL_DEF_(return,inline,set)
		saturate_f(auto &&x, auto &&...oo)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (provision::math::saturation_q<S_>) {
				return S_::template saturate_t<N_ism, N_car>::
					template method_f<Ns...>(XTAL_REF_(x), XTAL_REF_(oo)...);
			}
			XTAL_0IF (N_car == -0) {
				return XTAL_REF_(x);
			}
			XTAL_0IF (N_car == -1) {
				return XTAL_ALL_(x) {one};
			}
		}

	public:// CONSTRUCT
		using S_::S_;

	public:// TYPE

		using shape_parameter = typename metatype::shape_parameter;
		using   state_type = typename metatype::  state_type;
		using   order_attribute = typename metatype::  order_attribute;

	public:// OPERATE

		template <int N_ord=0, auto ...Ns> requires (1 <= N_ord)
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> s_gain
		,	atom::couple_q<unit_type[N_ord + 1]> auto &&scalars
		,	auto &&...oo
		)	const
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x);
			using U = unstruct_t<X>;
			using X_fit    = bond::fit<X>;
			using X_state_ = atom::couple_t<X[N_ord]>;
			using X_slope_ = atom::couple_t<X[N_ord]>;

			atom::math::dot_t<X[N_ord + 1]> outputs{};

			auto         outputs_ = outputs.self(cardinal_constant_t<N_ord>{});
			auto const   scalars_ = scalars.self(cardinal_constant_t<N_ord>{});
			auto const descalars_ = scalars_ - one;

			auto mem = S_::template memory<X_state_, X_slope_>();
			auto &states_ = get<0>(mem);
			auto &slopes_ = get<1>(mem);

		//	Initialize `slopes*`:
			slopes_.template blanket<1>();
			slopes_ *= scalars_;

		//	Initialize `outputs*`:
			get<0>(outputs) = get<0>(slopes_);
			bond::seek_until_f<N_ord - 1>([&] (auto I) XTAL_0FN {
				get<I + 1>(outputs_) = term_f(get<I + 1>(slopes_), get<I>(outputs_), s_gain);
			});
			get<N_ord>(outputs) = root_f<-1, (1)>(term_f(one, get<N_ord - 1>(outputs_), s_gain));

			int constexpr K_lim = provision::math::saturation_q<S_>;

		//	Integrate `states*` with `scalars*` and `outputs*`:
			bond::seek_until_f<K_lim + 1>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {get<N_ord>(outputs) *= x - dot_f(outputs_, states_);}
				XTAL_0IF (1 <= K) {get<N_ord>(outputs)  = x - dot_f(outputs_, slopes_);}

				bond::seek_until_f<-N_ord>([&] (auto I) XTAL_0FN {
					auto const o = get<I>(outputs_) = term_f(get<I>(states_), get<I + 1>(outputs), s_gain);
					if constexpr (K < K_lim) {
						get<I>(slopes_) = get<I>(descalars_) + saturate_f<-1, -1, Ns...>(o, oo...);
					}
				});
			});

		//	Finalize states and `outputs*`/`scalars*`:
			bond::seek_until_f<N_ord>([&] (auto I) XTAL_0FN {
				get<I>(states_) = term_f(-get<I>(states_), two, get<I>(outputs_));
			});
			static_assert(XTAL_ALL_(slopes_)::size() == XTAL_ALL_(outputs_)::size());
			slopes_ /= scalars_;
			reinterpret_cast<X_slope_ &>(outputs_) *= slopes_;//TODO: Make optional?

			return outputs.twin(constant_t<2>{});
		}
		template <int N_ord=0, auto ...Ns> requires (0 == N_ord)
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> s_gain
		,	atom::couple_q<unit_type[N_ord + 1]> auto &&scalars
		,	auto &&...oo
		)
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x);
			using U = unstruct_t<X>;
			return atom::math::dot_t<X[2]>{XTAL_REF_(x)};
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> s_gain
		,	atom::couple_q<unit_type[N_ord + 1]> auto &&scalars
		)
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x);
			using U = unstruct_t<X>;
			return method<N_ord, Ns...>(x, s_gain, XTAL_REF_(scalars), U{one});
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> s_gain
		,	unstruct_t<decltype(x)> s_damp
		,	auto &&...oo
		)
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x);
			using U = unstruct_t<X>;
			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return atom::math::dot_t<X[2]>{XTAL_REF_(x)};
			}
			XTAL_0IF (1 <= N_ord) {
				return method<N_ord, Ns...>(XTAL_REF_(x), s_gain
				,	[=] () XTAL_0FN -> atom::couple_t<U[N_ord + 1]> {
						auto const  &u = limit_f<[] XTAL_1FN_(to) (bond::fit<X>::minilon_f(1))>(s_damp);
						auto const u02 =             two*u , u04 = u02*two;
						auto const u12 = term_f(one, two,u), w24 = u02*u12;
						XTAL_IF0
						XTAL_0IF (1 == N_ord) {return {one, one};}
						XTAL_0IF (2 == N_ord) {return {one, u02, one};}
						XTAL_0IF (3 == N_ord) {return {one, u12, u12, one};}
						XTAL_0IF (4 == N_ord) {return {one, u04, w24, u04, one};}
					}()
				,	XTAL_REF_(oo)...
				);
			}
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> s_gain
		,	unstruct_t<decltype(x)> s_damp
		,	unstruct_t<decltype(x)> y_fade
		,	auto &&...oo
		)
		noexcept -> auto
		{
			using X =  XTAL_ALL_(x);
			using U = unstruct_t<X>;
			/**/
			return method<N_ord, Ns...>(XTAL_REF_(x), s_gain, s_damp, XTAL_REF_(oo)...)*
				atom::math::dot_t<U[2]>{term_f<-1, 2>(one, y_fade), y_fade};
			/*/
			return dot_f(method<N_ord, Ns...>(XTAL_REF_(x), s_gain, s_damp)
			,	W2{term_f<-1, 2>(one, y_fade), y_fade}
			);//TODO: Replace with `fade`.
			/***/
		}

		/*!
		\brief   Produces the `gain` and `damp` parameters from
		         the supplied `phason` and the `complex_field_q` `s`.
		
		Requires `1 <= Abs@s && Re@s <= 0 && 0 <= Im@s`.
		*/
		template <auto ...Ns>
		XTAL_DEF_(return,let)
		method(auto &&x
		,	atom::math::phason_simplex_q auto const &t_
		,	complex_field_q auto const &s
		,	auto &&...oo
		)
		noexcept -> auto
		{
			auto const &[s_re,  s_im] = destruct_f(XTAL_REF_(s));
			auto const  [s_a1, _s_a1] = roots_f<2>(square_f(s_re, s_im));

			auto const s_damp = s_im*_s_a1;
			auto       s_gain = t_(1)*s_a1;
			if constexpr (provision::math::prewarping_q<T_>) {
				s_gain *= pade::tangy_f<1, -1>(s_gain);
			}
			return method<Ns...>(XTAL_REF_(x), s_gain, s_damp, XTAL_REF_(oo)...);
		}
		/*!
		\brief   Produces the `gain` and `damp` parameters from
		         the supplied `phason` and `quason_q` triple `{beta, zeta, omega}`.
		\todo    Use the `phason` to reset the filter on discontinuity?
		*/		
		template <auto ...Ns>
		XTAL_DEF_(return,let)
		method(auto &&x
		,	atom::math::phason_simplex_q auto const &t_
		,	atom::math::quason_q<null_type[3]> auto const &o
		,	auto &&...oo
		)
		noexcept -> auto
		{
			auto const &[s_bend, s_damp, tau_abs] = o;
			auto s_gain = t_(1)/tau_abs;
			if constexpr (provision::math::prewarping_q<T_>) {
				s_gain *= pade::tangy_f<1, -1>(s_gain);
			}
			return method<Ns...>(XTAL_REF_(x), s_gain, s_damp, s_bend, XTAL_REF_(oo)...);
		}
		/*!
		\brief   Produces the `gain` and `damp` parameters from
		         the supplied `quason_q` triple `{beta, zeta, omega, phi}`.
		*/
		template <auto ...Ns>
		XTAL_DEF_(return,let)
		method(auto &&x
		,	atom::math::quason_q<null_type[4]> auto const &o
		,	auto &&...oo
		)
		noexcept -> auto
		{
			auto t_ = destruct_f<-1>(o);
			return method(XTAL_REF_(x), XTAL_MOV_(t_), XTAL_REF_(o).self(constant_t<-1>{}), XTAL_REF_(oo)...);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct context<process::math::zavalishin::filter<_s...>>
{
	using superkind = context<process::math::zavalishin::scaffold<_s...>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
	
	public:
		using S_::S_;

		template <extent_type N_mask=1>
		struct dispatch : bond::compose<void
		,	provision::voiced<void
			,	typename T_::   order_attribute::template dispatch<N_mask>
			>
		,	typename S_::template dispatch<N_mask>
		>
		{};
		template <extent_type N_mask=1>
		struct   attach : bond::compose<void
		,	provision::voiced<void
			,	typename T_::   stage_type::template  inspect<N_mask>
			,	typename T_::resample_type::template   attach<N_mask>
			>
		,	typename S_::template   attach<N_mask>
		>
		{};

	};
};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
