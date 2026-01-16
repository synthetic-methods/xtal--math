#pragma once
#include "./any.hh"

#include "./meta.hh"
#include "../../provision/prewarping.hh"
#include "../../provision/saturation.hh"
#include "../../atom/fourier/series.hh"

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
polymorphic and can accomodate up to the given cache-width (determined by `U_pole[N_pole][2]`).
For example, with base-types of `double` and `std::complex<double>` respectively,
the storage required is `16` and `32` bytes-per-pole.
*/
template <class ..._s>	struct  filter;
template <class ..._s>	concept filter_q = bond::tag_in_p<filter, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct filter
{
	using cotype = occur::context_t<filter>;

	using state_type = typename cotype::state_type;
	using slope_type = typename cotype::slope_type;

	using superkind = bond::compose<bond::tag<filter>
	,	provision::memorized<state_type, slope_type>
	,	meta<As...>
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

		using state_type      = typename cotype::state_type;
		using shape_parameter = typename cotype::shape_parameter;
		using order_attribute = typename cotype::order_attribute;

	private:
		XTAL_DEF_(set) M_ord = 0 + state_type::size();
		XTAL_DEF_(set) M_lim = 1 + state_type::size();

	public:// OPERATE

		template <int N_ord=0, auto ...Ns> requires (0 == N_ord)
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> s_gain
		,	atom::couple_q auto &&scalars_
	//	,	atom::couple_q<null_type[N_ord + 0]> auto &&scalars_
		,	auto &&...oo
		)
		noexcept -> atom::couple_t<XTAL_ALL_(x)[M_lim]>
		requires in_v<N_ord + 0, XTAL_ALL_(scalars_)::size()>
		{
			using X = XTAL_ALL_(x);
			return {XTAL_REF_(x)};
		}
		template <int N_ord=0, auto ...Ns> requires (1 <= N_ord)
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> s_gain
		,	atom::couple_q auto &&scalars_
	//	,	atom::couple_q<null_type[N_ord + 0]> auto &&scalars_
		,	auto &&...oo
		)	const
		noexcept -> atom::couple_t<XTAL_ALL_(x)[M_lim]>
		requires in_v<N_ord + 0, XTAL_ALL_(scalars_)::size()>
		{
			using X         = XTAL_ALL_(x);
			using X_fit     = bond::fit<X>;
			using X_state_  = atom::couple_t<X[N_ord]>;
			using X_slope_  = atom::couple_t<X[N_ord]>;

			atom::couple_t<X[M_lim]> outputs{};
			auto          outputs_ = outputs.self(cardinal_constant_t<N_ord>{});
		//	auto const    scalars_ = scalars.self(cardinal_constant_t<N_ord>{});

			auto mem = S_::template memory<X_state_, X_slope_>();
			auto &states_ = get<0>(mem);
			auto &slopes_ = get<1>(mem);

		//	Initialize `slopes*`:
			slopes_.template blanket<1>();
			slopes_ *= scalars_;

		//	Initialize `outputs*`:
			get<0>(outputs) = get<0>(slopes_);
			bond::seek_until_f<N_ord - 1>([&] (auto const I) XTAL_0FN {
				get<I + 1>(outputs_) = term_f(get<I + 1>(slopes_), get<I>(outputs_), s_gain);
			});
			get<N_ord>(outputs) = root_f<-1>(term_f(one,  get<N_ord - 1>(outputs_), s_gain));

			int constexpr K_lim = provision::math::saturation_q<S_>;

		//	Integrate `states*` with `scalars*` and `outputs*`:
			bond::seek_until_f<K_lim + 1>([&] (auto const K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {get<N_ord>(outputs) *= saturate_f<+1,-0>(x - dot_f(outputs_, states_), oo...);}
				XTAL_0IF (1 <= K) {get<N_ord>(outputs)  = saturate_f<+1,-0>(x - dot_f(outputs_, slopes_), oo...);}

				bond::seek_until_f<-N_ord>([&] (auto const I) XTAL_0FN {
					auto const o = get<I>(outputs_) = term_f(get<I>(states_), get<I + 1>(outputs), s_gain);
					if constexpr (K < K_lim) {
						get<I>(slopes_) = get<I>(scalars_) + saturate_f<-1,-1, Ns...>(o, oo...) - one;
					}
				});
			});

		//	Finalize states and `outputs*`/`scalars*`:
			bond::seek_until_f<N_ord>([&] (auto I) XTAL_0FN {
				get<I>(states_) = term_f(-get<I>(states_), two, get<I>(outputs_));
			});
			slopes_ /= scalars_;
			outputs_ *= slopes_;//TODO: Make optional?

			return outputs;
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> s_gain
		,	atom::couple_q auto &&scalars
		,	auto &&...oo
		)
		noexcept -> decltype(auto)
		requires in_v<N_ord + 1, XTAL_ALL_(scalars)::size()>
		{
			using X = XTAL_ALL_(x);
			using F = atom::math::fourier::series_t<X[N_ord]>;

			auto const [up, dn] = roots_f<N_ord>(get<N_ord>(scalars));
			auto rescalars = scalars.twin(constant_t<N_ord>{});
			rescalars *= F{dn};
			s_gain    *=  (up);
			return method<N_ord, Ns...>(x, s_gain, XTAL_MOV_(rescalars), XTAL_REF_(oo)...);
		}

	//	TODO: Move some of the specializations into `play`?

		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> s_gain
		,	atom::couple_q auto &&scalars
		)
		noexcept -> decltype(auto)
		{
			using X = XTAL_ALL_(x);
			return method<N_ord, Ns...>(x, s_gain, XTAL_REF_(scalars), unstruct_t<X>{one});
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> s_gain
		,	unstruct_t<decltype(x)> s_damp
		,	auto &&...oo
		)
		noexcept -> decltype(auto)
		{
			using X = XTAL_ALL_(x);
			return method<N_ord, Ns...>(XTAL_REF_(x), s_gain
			,	[=] () XTAL_0FN -> atom::couple_t<unstruct_t<X>[N_ord + 0]> {
					auto const u1_0 = limit_f<+1>(s_damp);
					auto const u2_0 = two*u1_0, u2_1 = u2_0 +  one;
					auto const u4_0 = two*u2_0, w4_2 = u2_0 * u2_1;
					XTAL_IF0
					XTAL_0IF (0 == N_ord) {return {/*one*/};}
					XTAL_0IF (1 == N_ord) {return {one/*, one*/};}
					XTAL_0IF (2 == N_ord) {return {one, u2_0/*, one*/};}
					XTAL_0IF (3 == N_ord) {return {one, u2_1, u2_1/*, one*/};}
					XTAL_0IF (4 == N_ord) {return {one, u4_0, w4_2, u4_0/*, one*/};}
				}	()
			,	XTAL_REF_(oo)...
			);
		}

		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> s_gain
		,	unstruct_t<decltype(x)> s_damp
		,	atom::math::dot_q auto &&y_dot
		,	auto &&...oo
		)
		noexcept -> objective_t<XTAL_ALL_(x)>
		{
			return XTAL_REF_(y_dot)*
				method<N_ord, Ns...>(XTAL_REF_(x), s_gain, s_damp, XTAL_REF_(oo)...);
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> s_gain
		,	unstruct_t<decltype(x)> s_damp
		,	unstruct_t<decltype(x)> y_fade
		,	auto &&...oo
		)
		noexcept -> objective_t<XTAL_ALL_(x)>
		{
			using X =  XTAL_ALL_(x);
			using U = unstruct_t<X>;
			return method<N_ord, Ns...>(XTAL_REF_(x), s_gain, s_damp,
				atom::math::dot_t<U[2]>{one - y_fade, y_fade}, XTAL_REF_(oo)...);
		}

		/*!
		\brief   Produces the `gain` and `damp` parameters from
		         the supplied `phason` and the `complex_field_q` `s`.
		
		Requires `1 <= Abs@s && Re@s <= 0 && 0 <= Im@s`.
		*/
		template <int N_ord=0, auto ...Ns>
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
	using superkind = context<process::math::zavalishin::meta<_s...>>;

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
