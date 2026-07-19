#pragma once
#include "./any.hh"

#include "./base.hh"
#include "../../atom/fourier/series.hh"
#include "../../scheme/zavalishin/distorted.hh"


XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Integrating filter using the Topology Preserving Transform (TPT).

Implementation outlined in _The Art of VA Synthesis_ by Vadim Zavalishin.

The non-linearity is supplied as a `process`-`template` using `scheme::math::zavalishin::distorted`.
The process must conform to the signature `<M_ism, M_car>`,
defining a stateless `method<N_var, ...>` within the `subtype`.

The parameters `M_ism` and `M_car` determine the type and return-value of shape, respectively.
`M_ism` is expected to yield convex/concave shapes for positive/negative valued,
and `M_car` is expected to return the slope when `== -1`.

\note    Despite the parameterization defined by `process::occurrence<filter<...>>`, `filter<...>::method` is
polymorphic and can accomodate up to the given cache-width (determined by `U_pole[N_pole][2]`).
For example, with base-types of `double` and `std::complex<double>` respectively,
the storage required is `16` and `32` bytes-per-pole.
*/
template <class ..._s>	struct  filter;
template <class ..._s>	concept filter_q = bond::tag_inner_p<filter, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct filter
{
	using    cotype = process::occurrence_t<filter>;
	using data_type = typename cotype::data_type;

	using superkind = bond::compose<bond::tag<filter>
	,	scheme::stashed<data_type[2]>
	,	base<_s...>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;

		template <int N_ism=0, int N_car=0, auto ...Ns>
		XTAL_DEF_(return,inline,set)
		f_distort(auto &&x, auto &&...oo)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (scheme::math::zavalishin::distorted_q<S_>) {
				using  distortion_type = typename S_::template distortion_t<N_ism, N_car>;
				return distortion_type{}.template method<Ns...>(XTAL_REF_(x), XTAL_REF_(oo)...);
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
		f_consort()
		noexcept -> int
		{
			return 0;
		}
		template <int N_ord, int N_ind, int N_dex=1> requires (1 <= N_ord)
		XTAL_DEF_(inline,set)
		f_consort()
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
		co_v = bond::operate_v<f_consort<Ns...>()>();


	public:// CONSTRUCT
		using S_::S_;

	public:// TYPE
		using data_type = typename S_::data_type;
		using order_attribute = typename S_::order_attribute;

	private:
		XTAL_DEF_(set) M_ord = 0 + data_type::size();
		XTAL_DEF_(set) M_lim = 1 + data_type::size();

	public:// FUSE
		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion> requires in_v<N_ion, -1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			using rewind_type = occur::rewind_t<XTAL_ALL_(o)>;
			using rewind_node = typename S_::template head_t<rewind_type>;
			XTAL_IF0
			XTAL_0IF (complete_q<rewind_node>) {
				auto const &q = self().template head<rewind_type>().wind();
				if (1 != q and q == o) {
					S_::stash(constant_t<>{});
				}
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}

	protected:// OPERATE
		using S_::self;

		template <int N_ord=M_ord, auto ...Ns> requires (0 == N_ord)
		XTAL_DEF_(inline,let)
		method_impl(auto &&x
		,	auto &&...
		)	const
		noexcept -> atom::couple_t<XTAL_ALL_(x)[M_lim]>
		{
			return {XTAL_REF_(x)};
		}
		template <int N_ord=M_ord, auto ...Ns> requires (1 <= N_ord)
		XTAL_DEF_(inline,let)
		method_impl(auto &&x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q<null_type[N_ord]> auto const &coeffs
	//	,	atom::couple_q<null_type[N_ord]> auto &&coeffs
		,	auto &&...oo
		)	const
		noexcept -> atom::couple_t<XTAL_ALL_(x)[M_lim]>
		{
			using  X        = XTAL_ALL_(x);
			using  X_fit    = bond::fit<X>;
			using  X_valued = atom::couple_t<X[M_lim]>;
			using  X_values = atom::couple_t<X[N_ord]>;
			using  X_slopes = atom::couple_t<X[N_ord]>;
			using  X_states = atom::couple_t<X[N_ord]>;
			auto constexpr U_cut_ = cut_e<[] XTAL_1FN_(to) (-X_fit::ratio_f(1, 3))>;
			auto constexpr X_cut_ = cut_e<[] XTAL_1FN_(to) (-X_fit::diplo_f(0x10))>;
			auto constexpr X_cut  = cut_f<[] XTAL_1FN_(to) (-X_fit::diplo_f(0x10))>;
		//	TODO: Using `X_cut_` fails under `clang/RELEASE`.

			X_valued valued;
			auto     values = valued.self(cardinal_constant_t<N_ord>{});// NOTE: `span`-based!
			auto     stash  = S_::template stash<X_states, X_slopes>();
			auto    &states = get<0>(stash);
			auto    &slopes = get<1>(stash);

			U_cut_(u_warp);
			u_warp *= pade::tangy_f< 1,-1>(u_warp);

		//	Denormalize `slopes` and initialize `valued`:
			bond::seek_to_e<N_ord>([&] (auto const I)
			XTAL_0FN -> void {
				auto &coeff = get<I>(coeffs);
				auto  slope = get<I>(slopes);
				slope *= coeff;
				slope += coeff;
				if constexpr (0 == I) get<I>(values) = term_f(slope);
				if constexpr (1 <= I) get<I>(values) = term_f(slope, u_warp, get<I - 1>(values));
			});
			get<N_ord>(valued) = root_f<-1>(term_f(one, u_warp, get<N_ord - 1>(valued)));

		//	Update `valued` by iterating `slopes`:
			auto constexpr  K_sat = scheme::math::zavalishin::distorted_q<S_>;
			auto constexpr  K_max = 2 * K_sat;
			auto constexpr  K_lim = 1 + K_max;
			bond::seek_to_e<K_lim>([&] (auto const K)
			XTAL_0FN -> void {
				auto &valve = get<N_ord>(valued);
				if constexpr (0 == K) valve *= X_cut(x - dot_f(values, states));
				if constexpr (1 <= K) valve  = X_cut(x - dot_f(values, slopes));
			//	X_cut_(valve);
				bond::seek_to_e<-N_ord>([&] (auto const I)
				XTAL_0FN -> void {
					auto &value = get<I>(valued) = X_cut(term_f(get<I>(states), get<I + 1>(valued), u_warp));
					auto &coeff = get<I>(coeffs);
					auto &slope = get<I>(slopes);
				//	X_cut_(value);
					if constexpr (0 == I) slope = XTAL_ALL_(slope){};
					if constexpr (1 <= I) slope = X_cut(f_distort<-1,-2, Ns...>(value, oo...));
					if constexpr (K <  K_max) {
						slope *= coeff;
						slope += coeff;
					}
				});
			});

		//	Accumulate `states` and carry distorted `slopes`:
			auto const &w_warp = u_warp *= two;
			bond::seek_to_e<-N_ord>([&] (auto const I)
			XTAL_0FN -> void {
				auto const &coeff = get<I + 0>(coeffs);
				auto const &slope = get<I + 0>(slopes);
				auto       &state = get<I + 0>(states);
				auto       &value = get<I + 0>(valued);
				auto       &valve = get<I + 1>(valued);
				state += valve*w_warp;
				if constexpr (1 <= I and 0 < K_max) {valve += coeff*slope*value;}
			});

		//	values *= coeffs;//TODO: Provide option?
			return valued;
		}

	public:// OPERATE
		template <int N_ord=M_ord, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	atom::math::dot_q auto &&o
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		{
			return XTAL_REF_(o)*
				method<N_ord, Ns...>(XTAL_REF_(x), XTAL_REF_(oo)...);
		}

		template <int N_ord=M_ord, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q<null_type[N_ord + 0]> auto &&coeffs_
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		requires in_v<N_ord + 0, XTAL_ALL_(coeffs_)::size()>
		{
			return method_impl<N_ord, Ns...>(XTAL_REF_(x),
				u_warp, XTAL_REF_(coeffs_), XTAL_REF_(oo)...);
		}

		template <int N_ord=M_ord, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q auto &&coeffs
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		requires in_v<N_ord + 1, XTAL_ALL_(coeffs)::size()>
		{
			using X = XTAL_ALL_(x);
			using F = atom::math::fourier::series_t<X[N_ord]>;

			auto const [up, dn] = roots_f<N_ord>(get<N_ord>(coeffs));
			auto recoeffs = coeffs.twin(cardinal_constant_t<N_ord>{});
			recoeffs *= F{dn};
			u_warp   *=  (up);
			return method<N_ord, Ns...>(x, u_warp, XTAL_MOV_(recoeffs), XTAL_REF_(oo)...);
		}

	//	TODO: Move some of the specializations into `play`?

		template <int N_ord=M_ord, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto const &x
		,	unstruct_t<decltype(x)> u_warp
		,	atom::couple_q auto &&coeffs
		)	const
		noexcept -> decltype(auto)
		{
			using X = XTAL_ALL_(x);
			return method<N_ord, Ns...>(x, u_warp, XTAL_REF_(coeffs), unstruct_t<X>{one});
		}
		template <int N_ord=M_ord, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x
		,	unstruct_t<decltype(x)> u_warp
		,	unstruct_t<decltype(x)> u_damp
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		{
			using X = XTAL_ALL_(x);
			using U = unstruct_t<X>;
			return method_impl<N_ord, Ns...>(XTAL_REF_(x), u_warp
			,	[a=cut_f<+1>(u_damp)]<auto ...I> (bond::seek_in_t<I...>)
					XTAL_0FN -> atom::couple_t<U[N_ord]> {
						return {term_f(co_v<N_ord, I, 0>, co_v<N_ord, I, 1>, a)...};
					}
				(bond::seek_to_t<N_ord>{})
			,	XTAL_REF_(oo)...
			);
		}

		/*!
		\brief   Produces the `gain` and `damp` parameters from
		         the supplied `phason` and the `complex_field_q` `s`.
		
		Requires `1 <= Abs@s && Re@s <= 0 && 0 <= Im@s`.
		*/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&x, atom::math::phason_q auto &&t_
		,	atom::quantity_multiplies_q auto &&s_
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		requires un_v<atom::math::phason_q<decltype(x)>>
		{
		//	auto &&[s_sig, s_mag] = XTAL_REF_(s_);
			auto const u_damp = std::imag(s_.signum());
			auto const u_warp = std::real(XTAL_REF_(t_)(1))*s_.template magnum<1>();
			return method<Ns...>(XTAL_REF_(x), u_warp, u_damp, XTAL_REF_(oo)...);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(atom::math::phason_q auto &&t_, auto &&x
		,	atom::quantity_multiplies_q auto &&s_
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		requires un_v<atom::math::phason_q<decltype(x)>>
		{
			return method<Ns...>(XTAL_REF_(x), XTAL_REF_(t_), XTAL_REF_(s_), XTAL_REF_(oo)...);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::process
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct occurrence<process::math::zavalishin::filter<_s...>>
{
	using superkind = occurrence<process::math::zavalishin::base<_s...>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
		using U_ = typename S_::data_type;
		using V_ = unstruct_t<U_>;
	
	public:
		using S_::S_;

		using  gain_parameter  = occur::inferred_t<_s..., union GAIN, V_>;
		using  damp_parameter  = occur::inferred_t<_s..., union DAMP, V_>;

		template <extent_type N_mask=1>
		struct dispatch : bond::compose<void
		,	typename T_::order_attribute::template dispatch<N_mask>
		,	typename S_::                 template dispatch<N_mask>
		>	{};

	};
};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
