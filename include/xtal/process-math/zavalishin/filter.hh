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
	using     zoom_type = occur::inferred_t<struct     ZOOM, typename bond::operating::alpha_type>;
	using   select_type = occur::inferred_t<struct   SELECT, bond::seek_t<0, 1>>;
	using    order_type = occur::inferred_t<struct    ORDER, bond::seek_t<0, 1, 2, 3>>;
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
			XTAL_0IF (is_q<    order_type, O>)             {S_::cache(0);}
			XTAL_0IF (is_q< topology_type, O>)             {S_::cache(0);}
			XTAL_0IF (occur::stage_q<      O>) {if (o == 0) S_::cache(0);}
			return S_::infuse(XTAL_REF_(o));
		}

		template <auto ...Ns>
		XTAL_DEF_(short)
		XTAL_LET method( auto &&x_input
		,	real_number_q auto   f_scale
		,	real_number_q auto   s_coeff
		)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(x_input)>;
			return method<Ns...>(XTAL_REF_(x_input), f_scale, s_coeff, _op::alpha_0);
		}
		template <int N_sel=0, int N_ord=0, int N_top=0, auto ...Ns>
		XTAL_DEF_(inline)
		XTAL_LET method( auto const &x_input
		,	real_number_q auto        f_scale
		,	real_number_q auto        s_coeff
		,	real_number_q auto        y_mix
		)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(x_input)>;

			XTAL_LET N  = N_ord + 1;
			XTAL_LET N_ = N_ord + 0;
			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return x_input;
			}
			XTAL_0IF (0 == N_top) {
				using U_input  = XTAL_ALL_(x_input);
				using U_exput  = XTAL_ALL_(x_input);
				using U_coeff  = typename _op::alpha_type;
				using W_exputs = algebra::sector_t<U_input[N ]>;
				using W_coeffs = algebra::sector_t<U_coeff[N ]>;
				
				union {W_exputs exputs; W_coeffs coeffs;} io{[=] () XTAL_0FN -> W_coeffs {
					auto constexpr _0 =  _op::alpha_0;
					auto constexpr _1 =  _op::alpha_1;
					auto constexpr _2 =  _op::alpha_2;
					auto const    a02 = term_f(_0, _2, s_coeff);
					auto const    a12 = term_f(_1, _2, s_coeff);
					XTAL_IF0
					XTAL_0IF (1 == N_ord) {return {_1, _1};}
					XTAL_0IF (2 == N_ord) {return {_1, a02, _1};}
					XTAL_0IF (3 == N_ord) {return {_1, a12, a12, _1};}
				}()};

				(void) edit<N_ord, N_top, Ns...>(io, x_input, f_scale);

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
		XTAL_LET edit(auto &io, auto const &x_input, real_number_q auto const &f_scale)
		noexcept -> void
		{
			using _op = bond::operate<decltype(x_input)>;
			XTAL_LET _1 = _op::alpha_1;
			XTAL_LET _2 = _op::alpha_2;
			XTAL_LET N_ =    N_ord + 0;
			XTAL_LET M_ =    N_ord - 1;

			using U_input   = XTAL_ALL_(x_input);
			using U_exput   = XTAL_ALL_(x_input);
			using U_coeff   = typename _op::alpha_type;
			using W_inputs_ = algebra::sector_t<U_input[N_]>;
			using W_exputs_ = algebra::sector_t<U_exput[N_]>;
			using W_coeffs_ = algebra::sector_t<U_coeff[N_]>;
			auto    cachet  = S_::template cache<W_inputs_, W_inputs_>();
			auto    coeffs  = io.coeffs; auto &coeffs_ = reinterpret_cast<W_coeffs_ &>(coeffs);
			auto   &exputs  = io.exputs; auto &exputs_ = reinterpret_cast<W_exputs_ &>(exputs);
			auto   &states_ = get<0>(cachet);
			auto   &slopes_ = get<1>(cachet);
		//	auto   [states_, slopes_] = S_::template cache<W_exputs_, W_exputs_>();//NOTE: Can't access from lambda...

			auto &sl_0 = get<0>(slopes_), sl_N = _1;
			auto &ex_0 = get<0>(exputs), &ex_N = get<N_>(exputs);

		//	Initialize coeffs `u*` and voltages `exputs*`:
		//	slopes_ *= coeffs_;
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				auto const &co_I = get<I>(coeffs_);
				auto       &sl_I = get<I>(slopes_);
				sl_I += sl_I == XTAL_ALL_(sl_I){0};// Reinitialize slopes to `1`...
				sl_I *= co_I;
			});

			ex_0 = sl_0;
			bond::seek_forward_f<M_>([&] (auto I) XTAL_0FN {
				get<I + 1>(exputs_) = term_f(get<I + 1>(slopes_), get<I>(exputs_), f_scale);
			});
			ex_N = root_f<-1, (4)>(term_f(sl_N, get<M_>(exputs_), f_scale));

			XTAL_LET I_lim = (size_type) N_lim != 0;// Only iterate when `sigmoid_t` is non-linear...

		//	Integrate states `s*` with coeffs/voltages:
			bond::seek_forward_f<1 + I_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {ex_N *= x_input - exputs_.product(states_);}
				XTAL_0IF (1 <= K) {ex_N  = x_input - exputs_.product(slopes_);}
				bond::seek_backward_f<N_>([&] (auto I) XTAL_0FN {
					auto const      &co_I = get<I>(coeffs_);
					auto const      &st_I = get<I>(states_);
					auto const       ex_I = term_f(st_I, get<I + 1>(exputs), f_scale);
					auto const       sl_I = sigmoid_t<-1, -1>::template function<N_lim, Ns...>(ex_I, co_I);//, get<I>(z_));
					get<I>(exputs_) = ex_I;
					get<I>(slopes_) = sl_I;
				});
			});

		//	Finalize states and voltages/coeffs:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				get<I>(states_) = term_f(-get<I>(states_), get<I>(exputs_), _2);
			});
			slopes_ /= coeffs_;
		//	exputs_ *= slopes_;// TODO: Make this line optional?
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
