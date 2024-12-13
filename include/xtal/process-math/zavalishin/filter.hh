#pragma once
#include "./any.hh"

#include "./sigmoid.hh"
#include "../root.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As> struct   filter;
template <typename ...As> using    filter_t = process::confined_t<filter<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct filter
{
	using     zoom_type = occur::inferred_t<struct     ZOOM, float>;
	using    order_type = occur::inferred_t<struct    ORDER, bond::seek_t<0, 1, 2, 3>>;
	using   select_type = occur::inferred_t<struct   SELECT, bond::seek_t<0, 1>>;
	using topology_type = occur::inferred_t<struct TOPOLOGY, bond::seek_t<0, 1>>;
	using    limit_type = occur::inferred_t<struct    LIMIT, bond::seek_t<0, 1>>;

	XTAL_SET N_order = order_type::cardinality() - 1;
	XTAL_SET N_cache = N_order<<1;
	using    U_cache = typename bond::operating::aphex_type;

//	TODO: Implement `maximum`? \

//	NOTE: Assuming a base-type of `std::complex`, \
	the minimum storage required is `(2)*(2*8) == 32` bytes-per-pole. \

	//\
	using subject = provision::context<void
	using subject = bond::compose<void
	,	typename     zoom_type::template   attach<>
	,	typename   select_type::template dispatch<>
	,	typename    order_type::template dispatch<>
	,	typename topology_type::template dispatch<>
	,	typename    limit_type::template dispatch<>// TODO: Replace with `sigmoid`'s parameters?
	>;
	using superkind = bond::compose<As...
	,	provision::cached<U_cache[N_cache]>
	,	provision::example<>
	,	subject
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <class O>
		XTAL_DEF_(short)
		XTAL_LET infuse(O &&o)
		noexcept -> signed
		{
			XTAL_IF0
			XTAL_0IF (is_q<O, order_type>) {
				S_::cache(0);
			}
			return S_::infuse(XTAL_REF_(o));
		}

		template <auto ...Ns>
		XTAL_DEF_(inline)
		XTAL_LET method(auto &&input
		,	real_number_q auto &&scale
		,	real_number_q auto &&coefficient
		)
		noexcept -> XTAL_ALL_(input)
		{
			using _op = bond::operate<decltype(input), decltype(scale)>;

			return method<Ns...>(XTAL_REF_(input), XTAL_REF_(scale), XTAL_REF_(coefficient), _op::alpha_0);
		}
		template <int N_sel=0, int N_ord=0, int N_top=0, auto ...Ns>
		XTAL_DEF_(inline)
		XTAL_LET method(auto &&input
		,	real_number_q auto &&scale
		,	real_number_q auto &&coefficient
		,	real_number_q auto &&mix
		)
		noexcept -> XTAL_ALL_(input)
		{
			using _op = bond::operate<decltype(input), decltype(scale)>;

			XTAL_LET N  = N_ord + 1;
			XTAL_LET N_ = N_ord + 0;
			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return XTAL_REF_(input);
			}
			XTAL_0IF (0 == N_top) {
				using U_input = XTAL_ALL_(input);
				using U_exput = XTAL_ALL_(input);
				using U_coeff = typename _op::alpha_type;

				using W_exput  = algebra::sector_t<U_input[N ]>;
				using W_exput_ = algebra::sector_t<U_input[N_]>;
				using W_coeff  = algebra::sector_t<U_coeff[N ]>;
				using W_coeff_ = algebra::sector_t<U_coeff[N_]>;
				
				union {W_exput exput; W_coeff coeff;} io{[=] () XTAL_0FN -> W_coeff {
					auto constexpr _0 =  _op::alpha_0;
					auto constexpr _1 =  _op::alpha_1;
					auto constexpr _2 =  _op::alpha_2;
					auto const    a02 = term_f(_0, _2, coefficient);
					auto const    a12 = term_f(_1, _2, coefficient);
					XTAL_IF0
					XTAL_0IF (1 == N_ord) {return {_1, _1};}
					XTAL_0IF (2 == N_ord) {return {_1, a02, _1};}
					XTAL_0IF (3 == N_ord) {return {_1, a12, a12, _1};}
				}()};

				(void) edit<N_ord, N_top, Ns...>(io, XTAL_REF_(input), XTAL_REF_(scale));

				XTAL_LET I_ =  static_cast<unsigned>(N_sel);
				XTAL_LET I0 = _std::countr_one(I_ >>  0) +  0, J0 = I0 + 1;
				XTAL_LET I1 = _std::countr_one(I_ >> J0) + J0, J1 = I1 + 1;
				return term_f(get<0>(io.exput), get<1>(io.exput), mix);
			}
			XTAL_0IF (1 == N_top) {
				return XTAL_REF_(input);
			}
		}
		template <int N_ord=0, int N_top=0, int N_lim=0, auto ...Ns> requires (1 <= N_ord and N_top == 0)
		XTAL_DEF_(inline)
		XTAL_LET edit(auto &io, auto &&input, real_number_q auto const &scale)
		noexcept
		{
			using _op = bond::operate<decltype(input), decltype(scale)>;
			XTAL_LET _1 = _op::alpha_1;
			XTAL_LET _2 = _op::alpha_2;
			XTAL_LET N_ =    N_ord + 0;
			XTAL_LET M_ =    N_ord - 1;

			using U_input = XTAL_ALL_(input);
			using U_exput = XTAL_ALL_(input);
			using U_coeff = typename _op::alpha_type;

			using W_input_ = algebra::sector_t<U_input[N_]>;
			using W_exput_ = algebra::sector_t<U_exput[N_]>;
			using W_coeff_ = algebra::sector_t<U_coeff[N_]>;
			auto    cache  = S_::template cache<W_input_, W_input_>();
			auto    coeff  = io.coeff; auto &coeff_ = reinterpret_cast<W_coeff_ &>(coeff);
			auto   &exput  = io.exput; auto &exput_ = reinterpret_cast<W_exput_ &>(exput);
			auto   &slope_ = get<0>(cache);
			auto   &state_ = get<1>(cache);
		//	auto   [state_, slope_] = S_::template cache<W_exput_, W_exput_>();//NOTE: Can't access from lambda...

			auto &sl_0 = get<0>(slope_), sl_N = _1;
			auto &ex_0 = get<0>(exput), &ex_N = get<N_>(exput);

		//	Initialize coeff `u*` and voltages `exput*`:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				auto const &co_I = get<I>(coeff_);
				auto       &sl_I = get<I>(slope_);
				sl_I += sl_I == U_input{};// Reinitialize slopes to `1`...
				sl_I *= co_I;
			});
		//	slope_ *= coeff_;

			ex_0 = sl_0;
			bond::seek_forward_f<M_>([&] (auto I) XTAL_0FN {
				get<I + 1>(exput_) = term_f(get<I + 1>(slope_), get<I>(exput_), scale);
			});
			ex_N = root_f<-1, (4)>(term_f(sl_N, get<M_>(exput_), scale));

			XTAL_LET I_lim = (size_type) N_lim != 0;// Only iterate when `sigmoid_t` is non-linear...

		//	Integrate states `s*` with coeff/voltages:
			bond::seek_forward_f<1 + I_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {ex_N *= input - exput_.product(state_);}
				XTAL_0IF (1 <= K) {ex_N  = input - exput_.product(slope_);}
				bond::seek_backward_f<N_>([&] (auto I) XTAL_0FN {
					auto const      &co_I = get<I>(coeff_);
					auto const      &st_I = get<I>(state_);
					auto const       ex_I = term_f(st_I, get<I + 1>(exput), scale);
					auto const       sl_I = sigmoid_t<-1,-1>::template function<N_lim, Ns...>(ex_I, co_I);//, get<I>(z_));
					get<I>(exput_) = ex_I;
					get<I>(slope_) = sl_I;
				});
			});

		//	Finalize states and voltages/coeff:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				get<I>(state_) = term_f(-get<I>(state_), get<I>(exput_), _2);
			});
			slope_ /= coeff_;
		//	exput_ *= slope_;// TODO: Make this line optional?

			return exput;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
