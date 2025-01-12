#pragma once
#include "./any.hh"

#include "../../provision/shaper.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Integrating filter based on the Topology Preserving Transform (TPT), \
described in _The Art of VA Synthesis_ by Vadim Zavalishin. \

///\
The non-linearity is supplied as a `process`-`template` using `provision::shaper`. \
The process must conform to the signature `<M_ism, M_car>`, \
defining a stateless `function<N_var, ...>` within the `subtype`. \

///\
The parameters `M_ism` and `M_car` determine the type and return-value of curve, respectively. \
`M_ism` is expected to yield convex/concave curves for positive/negative values, \
and `M_car` is expected to return the slope when `== -1`. \

///\note\
The filter is preconfigured to support `std::complex<double>[4]` poles. \
With base-types of `double` and `std::complex<double>` respectively, \
the storage required is `16` and `32` bytes-per-pole. \

template <typename ...As> struct   filter;
template <typename ...As> using    filter_t = process::confined_t<filter<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class U_pole, int N_pole>
struct filter<U_pole[N_pole]>
{
	using     zoom_type = occur::inferred_t<struct     ZOOM, typename bond::operating::alpha_type>;
	using   select_type = occur::inferred_t<struct   SELECT, unsigned int, bond::word<2>>;
	using    order_type = occur::inferred_t<struct    ORDER, unsigned int, bond::seek_s<N_pole + 1>>;
	using topology_type = occur::inferred_t<struct TOPOLOGY, unsigned int, bond::word<2>>;
	using    limit_type = occur::inferred_t<struct    LIMIT, unsigned int, bond::word<2>>;

	using superkind = bond::compose<bond::tag<filter_t>
	,	provision::cached<U_pole[N_pole << 1]>
	,	provision::example<>
	,	typename     zoom_type::template   attach<>
	,	typename   select_type::template dispatch<>
	,	typename    order_type::template dispatch<>
	,	typename topology_type::template dispatch<>
	,	typename    limit_type::template dispatch<>// TODO: Replace with `shape`'s parameters?
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <class X>
		XTAL_DEF_(short)
		XTAL_LET infuse(X &&x)
		noexcept -> signed
		{
			XTAL_IF0
			XTAL_0IF (in_q<X, order_type, topology_type>) {
				S_::cache(constant_t<0>{});
			}
			XTAL_0IF (occur::stage_q<X>) {
				if (x == 0) S_::cache(constant_t<0>{});
			}
			return S_::infuse(XTAL_REF_(x));
		}

		template <auto ...Ns>
		XTAL_DEF_(short)
		XTAL_LET method( auto &&x_input
		,	real_variable_q auto   g_scale
		,	real_variable_q auto   s_coeff
		)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(x_input)>;
			return method<Ns...>(XTAL_REF_(x_input), g_scale, s_coeff, _op::alpha_0);
		}
		template <int N_sel=0, int N_ord=0, int N_top=0, auto ...Ns>
		XTAL_DEF_(inline)
		XTAL_LET method( auto const &x_input
		,	real_variable_q auto g_scale
		,	real_variable_q auto s_coeff
		,	real_variable_q auto y_mix
		)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(x_input)>;

			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return x_input;
			}
			XTAL_0IF (0 == N_top) {
				using U_input  = XTAL_ALL_(x_input);
				using U_exput  = XTAL_ALL_(x_input);
				using U_coeff  = typename _op::alpha_type;
				using U_exputs = arrange::couple_t<U_input[N_ord + 1]>;
				using U_coeffs = arrange::couple_t<U_coeff[N_ord + 1]>;
				
				union {U_exputs exputs; U_coeffs coeffs;} io{[=] ()
				XTAL_0FN -> U_coeffs {
					XTAL_LET K_1 = _op::alpha_f(1);
					XTAL_LET K_2 = _op::alpha_f(2);
					auto const  &u = s_coeff;
					auto const u02 =             K_2*u , u04 = u02*K_2;
					auto const u12 = term_f(K_1, K_2,u), w24 = u02*u12;
					XTAL_IF0
					XTAL_0IF (1 == N_ord) {return {K_1, K_1};}
					XTAL_0IF (2 == N_ord) {return {K_1, u02, K_1};}
					XTAL_0IF (3 == N_ord) {return {K_1, u12, u12, K_1};}
					XTAL_0IF (4 == N_ord) {return {K_1, u04, w24, u04, K_1};}
				}()};

				(void) edit<N_ord, N_top, Ns...>(io, x_input, g_scale);// io.coeffs -> io.exputs

				XTAL_LET I_ =  static_cast<unsigned>(N_sel);
				XTAL_LET I0 = _std::countr_one(I_ >>  0) +  0, J0 = I0 + 1;
				XTAL_LET I1 = _std::countr_one(I_ >> J0) + J0, J1 = I1 + 1;
				return term_f(get<0>(io.exputs), get<1>(io.exputs), y_mix);
			}
			XTAL_0IF (1 == N_top) {
				return x_input;
			}
		}
		template <int N_ord=0, int N_top=0, int N_lim=0, auto ...Ns> requires (1 <= N_ord and N_top == 0)
		XTAL_DEF_(inline)
		XTAL_LET edit(auto &io, auto const &x_input, real_variable_q auto const &g_scale)
		noexcept -> void
		{
			using _op = bond::operate<decltype(x_input)>;

			auto constexpr K_1 = _op::alpha_f(1);
			auto constexpr K_2 = _op::alpha_f(2);
			auto constexpr N_  =     N_ord - (0);
			auto constexpr M_  =     N_ord - (1);

			using U_input   = XTAL_ALL_(x_input);
			using U_exput   = XTAL_ALL_(x_input);
			using U_coeff   = typename _op::alpha_type;
			using U_inputs_ = arrange::couple_t<U_input[N_]>;
			using U_exputs_ = arrange::couple_t<U_exput[N_]>;
			using U_coeffs_ = arrange::couple_t<U_coeff[N_]>;

			auto    cachet  = S_::template cache<U_inputs_, U_inputs_>();
			auto    coeffs  = io.coeffs; auto &coeffs_ = reinterpret_cast<U_coeffs_ &>(coeffs);
			auto   &exputs  = io.exputs; auto &exputs_ = reinterpret_cast<U_exputs_ &>(exputs);
			auto   &slopes_ = get<0>(cachet);
			auto   &states_ = get<1>(cachet);
		//	auto   [slopes_, states_] = S_::template cache<U_exputs_, U_exputs_>();//NOTE: Can't access from lambda...

			auto   &sl_0    = get<0>(slopes_), sl_N = K_1;
			auto   &ex_0    = get<0>(exputs), &ex_N = get<N_>(exputs);

		//	Initialize `coeffs*`:
			slopes_.template unzero<1>();
			slopes_ *= coeffs_;

		//	Initialize `exputs*`:
			ex_0 = sl_0;
			bond::seek_forward_f<M_>([&] (auto I) XTAL_0FN {
				get<I + 1>(exputs_) = term_f(get<I + 1>(slopes_), get<I>(exputs_), g_scale);
			});
			ex_N = root_f<-1, (4)>(term_f(sl_N, get<M_>(exputs_), g_scale));

			XTAL_LET K_lim = provision::shaper_q<S_> and above_p<0, N_lim>;

		//	Integrate `states*` with `coeffs*` and `exputs*`:
			bond::seek_forward_f<1 + 1*K_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {ex_N *= x_input - exputs_.product(states_);}
				XTAL_0IF (1 <= K) {ex_N  = x_input - exputs_.product(slopes_);}

				bond::seek_backward_f<N_>([&] (auto I) XTAL_0FN {
					XTAL_IF0
					XTAL_0IF (0 == K_lim) {
						get<I>(exputs_)  = term_f(get<I>(states_), get<I + 1>(exputs), g_scale);
					}
					XTAL_0IF (1 == K_lim) {
						auto const exput = term_f(get<I>(states_), get<I + 1>(exputs), g_scale);
						auto const slope = S_::template shape_t<-1, -1>::template function<N_lim, Ns...>(exput) + (get<I>(coeffs_) - one);
						get<I>(exputs_) = exput;
						get<I>(slopes_) = slope;
					}
				});
			});

		//	Finalize states and `exputs*`/`coeffs*`:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				get<I>(states_) = term_f(-get<I>(states_), get<I>(exputs_), K_2);
			});
			slopes_ /= coeffs_;
		//	exputs_ *= slopes_;// TODO: Make this line optional?
		}

	};
};
template <>
struct filter<>
:	filter<typename bond::operating::aphex_type[4]>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
