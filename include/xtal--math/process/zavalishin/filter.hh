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
			XTAL_0IF (N_car == -0) {return XTAL_REF_(x)       ;}// f[x]
			XTAL_0IF (N_car == -1) {return XTAL_ALL_(x) {one };}// f[x]/x
			XTAL_0IF (N_car == -2) {return XTAL_ALL_(x) {zero};}// f[x]/x - 1
		}

		/*!
		Produces `CoefficientList[#, {s, r}][[1 + ind, 1 + dex]] &` in the expansion of
		```m
		(1 + s)                             (* 1 == ord *)
		(1 + 2 r s + s^2)*(1 + s)^(ord - 2) (* 2 <= ord *)
		```
		*/
		template <int N_ord, int N_ind, int N_dex=1> requires (0 == N_ord)
		XTAL_DEF_(inline,set)
		fact_f()
		noexcept -> int
		{
			return 0;
		}
		template <int N_ord, int N_ind, int N_dex=1> requires (1 <= N_ord)
		XTAL_DEF_(inline,set)
		fact_f()
		noexcept -> int
		{
			auto constexpr I_ind = modulo_v<N_ord, N_ind*sign_v<N_ord|1 - N_ind*2>>;
			auto constexpr I_ord = N_ord - 3;
			XTAL_IF0
			XTAL_0IF (0 == N_dex) {
				int o = 1;
				for (int i{I_ind}; ~--i;) {
					o *= I_ord - i;
					o += factorial_f(I_ind, i);
				}
				return o/factorial_f(I_ind);
			}
			XTAL_0IF (1 == N_dex and 0 == N_ind) {
				return 0;
			}
			XTAL_0IF (1 == N_dex and 1 <= N_ind) {
				int constexpr O = I_ord + 1;
				int constexpr I = I_ind - 1;
				static_assert(I <= O);
				return two*factorial_f(O, I)/factorial_f(O - I);
			}
		}
		template <int ...Ns>
		XTAL_DEF_(set)
		fact_v = bond::operate_v<fact_f<Ns...>()>();


	public:// CONSTRUCT
		using S_::S_;

	public:// TYPE

		using state_type      = typename S_::state_type;
		using shape_parameter = typename S_::shape_parameter;
		using order_attribute = typename S_::order_attribute;

	private:
		XTAL_DEF_(set) M_ord = 0 + state_type::size();
		XTAL_DEF_(set) M_lim = 1 + state_type::size();

	protected:// OPERATE

		template <int N_ord=0, auto ...Ns> requires (0 == N_ord)
		XTAL_DEF_(inline,let)
		methodology(auto const &x
		,	auto &&...
		)
		noexcept -> atom::couple_t<XTAL_ALL_(x)[M_lim]>
		{
			return {XTAL_REF_(x)};
		}
		template <int N_ord=0, auto ...Ns> requires (1 <= N_ord)
		XTAL_DEF_(inline,let)
		methodology(auto const &x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q auto const &coeffs_
	//	,	atom::couple_q<null_type[N_ord]> auto &&coeffs_
		,	auto &&...oo
		)	const
		noexcept -> atom::couple_t<XTAL_ALL_(x)[M_lim]>
		requires in_v<N_ord, XTAL_ALL_(coeffs_)::size()>
		{
			using  X         = XTAL_ALL_(x);
			using  U         = unstruct_t<X>;
			using  X_values  = atom::couple_t<X[M_lim]>;
			using  X_slopes_ = atom::couple_t<X[N_ord]>;
			using  X_states_ = atom::couple_t<X[N_ord]>;

			X_values values;
			auto     values_ = values.self(cardinal_constant_t<N_ord>{});// NOTE: `span`-based!
			auto     mem     = S_::template memory<X_states_, X_slopes_>();
			auto    &states_ = get<0>(mem);
			auto    &slopes_ = get<1>(mem);

		//	Denormalize `slopes` and initialize `values`:
			bond::seek_until_f<N_ord>([&] (auto const I) XTAL_0FN {
				auto const &coeff = get<I>(coeffs_);
				auto        slope = get<I>(slopes_);
				//\
				slope *= fact_v<N_ord, I, 1>;
				slope *= coeff;
				slope += coeff;
				if constexpr (0 == I) get<I>(values_) = term_f(slope);
				if constexpr (1 <= I) get<I>(values_) = term_f(slope, u_warp, get<I - 1>(values_));
			});
			get<N_ord>(values) = root_f<-1>(term_f(one, u_warp, get<N_ord - 1>(values)));

		//	Update `values` by iterating `slopes`:
			auto constexpr K_sat = provision::math::saturation_q<S_>;
			auto constexpr K_max = term_f(0, 2, K_sat);
			auto constexpr K_lim = term_f(1, 1, K_max);
			bond::seek_until_f<K_lim>([&] (auto const K) XTAL_0FN {
				if constexpr (0 == K) {get<N_ord>(values) *= x - dot_f(values_, states_);}
				if constexpr (1 <= K) {get<N_ord>(values)  = x - dot_f(values_, slopes_);}
				bond::seek_until_f<-N_ord>([&] (auto const I) XTAL_0FN {
					auto const &value = get<I>(values_) = term_f(get<I>(states_), get<I + 1>(values), u_warp);
					auto const &coeff = get<I>(coeffs_);
					auto       &slope = get<I>(slopes_);
					if constexpr (0 <= I) {slope  = zero;}
					if constexpr (1 <= I) {slope  = saturate_f<-1,-2, Ns...>(value, oo...);}
					if constexpr (K <  K_max) {
						//\
						slope *= fact_v<N_ord, I, 1>;
						slope *= coeff;
						slope += coeff;
					}
				});
			});

		//	Accumulate `states` and carry saturated `slopes`:
			bond::seek_until_f<-N_ord>([&, w_warp=u_warp*two] (auto const I) XTAL_0FN {
				auto const &coeff  = get<I + 0>(coeffs_);
				auto const &slope  = get<I + 0>(slopes_);
				auto       &state  = get<I + 0>(states_);
				auto       &value  = get<I + 0>(values);
				auto       &value_ = get<I + 1>(values);
				state += w_warp*value_;
				if constexpr (1 <= I and 0 < K_max) {
					value_ += value*coeff*slope;
				}
			});

		//	values_ *= coeffs_;//TODO: Provide option?
			return values;
		}


	public:// OPERATE

		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q auto &&coeffs_
	//	,	atom::couple_q<null_type[N_ord + 0]> auto &&coeffs_
		,	auto &&...oo
		)
		noexcept -> atom::couple_t<XTAL_ALL_(x)[M_lim]>
		requires in_v<N_ord + 0, XTAL_ALL_(coeffs_)::size()>
		{
			return methodology<N_ord, Ns...>(XTAL_REF_(x), u_warp, XTAL_REF_(coeffs_), XTAL_REF_(oo)...);
		}


		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q auto &&coeffs
		,	auto &&...oo
		)
		noexcept -> decltype(auto)
		requires in_v<N_ord + 1, XTAL_ALL_(coeffs)::size()>
		{
			using X = XTAL_ALL_(x);
			using F = atom::math::fourier::series_t<X[N_ord]>;

			auto const [up, dn] = roots_f<N_ord>(get<N_ord>(coeffs));
			auto recoeffs = coeffs.twin(constant_t<N_ord>{});
			recoeffs *= F{dn};
			u_warp   *=  (up);
			return method<N_ord, Ns...>(x, u_warp, XTAL_MOV_(recoeffs), XTAL_REF_(oo)...);
		}

	//	TODO: Move some of the specializations into `play`?

		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q auto &&coeffs
		)
		noexcept -> decltype(auto)
		{
			using X = XTAL_ALL_(x);
			return method<N_ord, Ns...>(x, u_warp, XTAL_REF_(coeffs), unstruct_t<X>{one});
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> u_warp
		,	unstruct_t<decltype(x)> u_damp
		,	auto &&...oo
		)
		noexcept -> decltype(auto)
		{
			using X  = XTAL_ALL_(x);
			using U  = unstruct_t<X>;
			return methodology<N_ord, Ns...>(XTAL_REF_(x), u_warp
			,	[a=limit_f<+1>(u_damp)]<auto ...I> (bond::seek_t<I...>)
					XTAL_0FN -> atom::couple_t<U[N_ord]> {
						return {term_f(fact_v<N_ord, I, 0>, fact_v<N_ord, I, 1>, a)...};
					}
				(bond::seek_s<N_ord>{})
			,	XTAL_REF_(oo)...
			);
		}

		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> u_warp
		,	unstruct_t<decltype(x)> u_damp
		,	atom::math::dot_q auto &&y_dot
		,	auto &&...oo
		)
		noexcept -> objective_t<XTAL_ALL_(x)>
		{
			return XTAL_REF_(y_dot)*
				method<N_ord, Ns...>(XTAL_REF_(x), u_warp, u_damp, XTAL_REF_(oo)...);
		}
		template <int N_ord=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> u_warp
		,	unstruct_t<decltype(x)> u_damp
		,	unstruct_t<decltype(x)> y_fade
		,	auto &&...oo
		)
		noexcept -> objective_t<XTAL_ALL_(x)>
		{
			using X =  XTAL_ALL_(x);
			using U = unstruct_t<X>;
			auto const y_up =       y_fade, v_up = root_f<2>(y_up);
			auto const y_dn = one - y_fade, v_dn = root_f<2>(y_dn);
			return method<N_ord, Ns...>(XTAL_REF_(x), u_warp, u_damp,
				atom::math::dot_t<U[2]>{y_dn, y_up}, XTAL_REF_(oo)...);
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

			auto const u_damp = s_im*_s_a1;
			auto       u_warp = t_(1)*s_a1;
			if constexpr (provision::math::prewarping_q<T_>) {
				u_warp *= pade::tangy_f<1, -1>(u_warp);
			}
			return method<Ns...>(XTAL_REF_(x), u_warp, u_damp, XTAL_REF_(oo)...);
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
