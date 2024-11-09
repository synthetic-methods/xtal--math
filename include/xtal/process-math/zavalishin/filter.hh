#pragma once
#include "./any.hh"

#include "./sigmoid.hh"
#include "../root.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As> XTAL_TYP filter;
template <typename ...As> XTAL_USE filter_t = process::confined_t<filter<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct filter
:	process::link<filter<>, As...>
{
};
template <>
struct filter<>
{
	using  alpha_type = typename bond::operating::alpha_type;
	using  aphex_type = typename bond::operating::aphex_type;

	using   zoom_type = occur::inferred_t<XTAL_TYP   ZOOM, float>;
	using  order_type = occur::inferred_t<XTAL_TYP  ORDER, bond::seek_t<0, 1, 2, 3>>;
	using  limit_type = occur::inferred_t<XTAL_TYP  LIMIT, bond::seek_t<0, 1>>;
	using select_type = occur::inferred_t<XTAL_TYP SELECT, bond::seek_t<0, 1>>;

	XTAL_SET N_order = order_type::cardinality() - 1;
	XTAL_SET N_cache = N_order*2;

//	TODO: Implement `maximum`? \

//	NOTE: Assuming a base-type of `std::complex`, \
	the minimum storage required is `(2)*(2*8) == 32` bytes-per-pole. \

	//\
	using subject = resource::invoice<void
	using subject = bond::compose<void
	,	typename   zoom_type::template   attach<>
	,	typename  limit_type::template dispatch<>
	,	typename  order_type::template dispatch<>
	,	typename select_type::template dispatch<>
	>;
	using superkind = bond::compose<void
	,	resource::cached<aphex_type[N_cache]>
	,	resource::example<>
	,	subject
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		XTAL_TNX infuse(auto &&o)
		XTAL_0EX
		{
			XTAL_IF0
			XTAL_0IF (is_q<decltype(o), order_type>) {
			//	Reset? Damb...
			}
			return S_::infuse(XTAL_REF_(o));
		}
		template <int N_lim=0, int N_ord=0, int N_sel=0>
		XTAL_DEF_(inline)
		XTAL_LET method(auto &&x
		,	real_number_q auto &&o
		,	real_number_q auto &&a
		)
	//	XTAL_0EX -> algebra::scalar_t<XTAL_ALL_(x)[N_ord]>
		XTAL_0EX -> XTAL_ALL_(x)
		{
			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return XTAL_REF_(x);
			}
			XTAL_0IF (0 == N_sel) {
				using horner::term_f;
				auto constexpr N  = N_ord + 1;
				auto constexpr N_ = N_ord + 0;

				using X  = XTAL_ALL_(x); using _op = bond::operate<X>;
				using A  = algebra::scalar_t<typename _op::alpha_type[N ]>;
				using A_ = algebra::scalar_t<typename _op::alpha_type[N_]>;
				
				auto constexpr _1 = _op::alpha_1;
				auto constexpr _2 = _op::alpha_2;

				XTAL_IF0
				XTAL_0IF (1 == N_ord) {
					return method<N_lim, N_ord, N_sel>(XTAL_REF_(x), XTAL_REF_(o), A_{_1});
				}
				XTAL_0IF (2 == N_ord) {
					auto const a1 = _2*a;
					return method<N_lim, N_ord, N_sel>(XTAL_REF_(x), XTAL_REF_(o), A_{_1, a1});
				}
				XTAL_0IF (3 == N_ord) {
					auto const a1 = term_f(_1, _2, a);
					return method<N_lim, N_ord, N_sel>(XTAL_REF_(x), XTAL_REF_(o), A_{_1, a1, a1});
				}
			}
			XTAL_0IF (1 == N_sel) {
				return XTAL_REF_(x);
			}
		}
		template <int N_lim=0, int N_ord=0, int N_sel=0> requires (1 <= N_ord and N_sel == 0)
		XTAL_DEF_(inline)
		XTAL_LET method(auto &&x
		,	real_number_q auto const &o
		,	algebra::scalar_q auto const &a_
		)
		XTAL_0EX -> auto
		{
			using horner::term_f;
			auto constexpr N  = N_ord + 1;
			auto constexpr N_ = N_ord + 0;
			auto constexpr M_ = N_ord - 1;

			using X  = XTAL_ALL_(x); using _op = bond::operate<X>;
			using W  = algebra::scalar_t<X[N ]>;
			using W_ = algebra::scalar_t<X[N_]>;

		//	static_assert(is_q<W , decltype(a )>);
			static_assert(is_q<W_, decltype(a_)>);

		//	algebra::series_t<X[N_]> z_(self().template head<zoom_type>());
		//	auto [s_, u_] = S_::template cache<W_, W_>();//NOTE: Can't access from lambda...
			auto cache = S_::template cache<W_, W_>();
			W w;
			auto &w_ = reinterpret_cast<W_ &>(w);
			auto &u_ = get<0>(cache);
			auto &s_ = get<1>(cache);

			auto &u0 = get<0>(u_), uN = _op::alpha_1;
			auto &w0 = get<0>(w), &wN = get<N_>(w);

		//	Initialize coefficients `u*` and voltages `w*`:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				auto &uI = get<I>(u_);
				uI += uI == X{};// Force non-zero... (annoying)
				uI *= get<I>(a_);
			});
		//	u_ *= a_;
			w0  = u0;
			bond::seek_forward_f<M_>([&] (auto I) XTAL_0FN {
				get<I + 1>(w_) = term_f(get<I + 1>(u_), get<I>(w_), o);
			});
			wN = root_f<-1, (4)>(term_f(uN, get<M_>(w_), o));

		//	Integrate states `s*` with coefficients/voltages:
			bond::seek_forward_f<1 + N_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {wN *= x - w_.product(s_);}
				XTAL_0IF (1 <= K) {wN  = x - w_.product(u_);}
				bond::seek_backward_f<N_>([&] (auto I) XTAL_0FN {
					auto const  &aI = get<I>(a_);
					auto const  &sI = get<I>(s_);
					auto const   wI = term_f(sI, get<I + 1>(w), o);
					auto const   uI = sigmoid_t<-1,-1>::template function<N_lim>(wI, aI);//, get<I>(z_));
					get<I>(w_) = wI;
					get<I>(u_) = uI;
				});
			});
		//	Finalize states and voltages/coefficients:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				get<I>(s_) = term_f(-get<I>(s_), get<I>(w_), _op::alpha_2);
			});
		//	s_  = term_f(-s_, w_, _op::alpha_2);
			u_ /= a_;
		//	w_ *= u_;// TODO: Make this line optional?

		//	return y;
			return get<0>(w);
		}
		/*/
		template <int N_lim=0, int N_ord=0, int N_sel=0> requires (0 == N_ord and 0 == N_sel)
		XTAL_DEF_(return)
		XTAL_LET method(auto &&x, auto &&f)
		XTAL_0EX
		{
			using namespace horner;
			using V = XTAL_ALL_(x); using _op = bond::operate<V>;
			using A = typename bond::operate<V>::alpha_type;
			
			using W_ = algebra::scalar_t<V[3]>;
			using V_ = algebra::scalar_t<V[2]>;
			using A_ = algebra::scalar_t<A[2]>;
			
			auto [w_, a_] = S_::template cache<V_, A_>();
			A const &a0 = a_[0];
			A const &a1 = a_[1];
			A const  a2 =    1 ;

			A const w0 =        a0        ;// a_[0]
			A const w1 = term_f(a1, w0, f);// a_[1] + f %
			A const w2 = term_f(a2, w1, f);// a_[2] + f %
			W_ y_{V{}, V{}, term_f<-1>(term_f<-1>(x, w_[0], w0), w_[1], w1)/w2};

		//	solve<Is...>(y_, x, f);// TODO: Save current `y_[1]`!
			
			y_[1] = term_f(w_[1], y_[2], f); w_[1] = term_f(y_[1], y_[2], f);
			y_[0] = term_f(w_[0], y_[1], f); w_[0] = term_f(y_[0], y_[1], f);

			return y_;
		}
		/***/
		/*/
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&x, auto &&f, auto &&_R)
		XTAL_0EX
		{
			using V = XTAL_ALL_(x); using _op = bond::operate<V>;
			using A = typename _op::alpha_type;

			using V_ = algebra::scalar_t<V[2]>;
			using A_ = algebra::scalar_t<A[2]>;

			auto [w_, a_] = S_::template cache<V_, A_>();

		//	a1 = 2*cos(1.5707963267948966*alpha);
			new (&a_) A_{_op::alpha_1, _op::alpha_2*_R};
			return method<Is...>(XTAL_REF_(x), f);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&x, auto &&f, complex_number_q auto &&_S)
		XTAL_0EX
		{
			//\todo\
			Use `bicycle` to provide the base frequency, \
			combining relative-frequency and resonance as the complex `s`, \
			where frequency is given by `Abs[s]`, and `Q = Sqrt[# - ¼]@?`. \

			using V = XTAL_ALL_(x); using _op = bond::operate<V>;
			using A = typename _op::alpha_type;

			auto constexpr N_lnH = (A) -0.6931471805599453094172321214581765681e+0L;
			
			return method(XTAL_REF_(x), XTAL_REF_(f), exp(N_lnH*XTAL_REF_(_S).imag()));
		}
		/***/
		/*/
		template <auto ...Is>
		XTAL_DEF_(inline)
		XTAL_LET solve(auto &y_, auto const &x, auto const &f)
		XTAL_0EX -> void
		{
//		With:
//			f(y[1]) = y1*(a_[1] - 2) + 2*F[y1, o, v]
//		
//		Fixed-point:
//			y1 = term_f< 1>(w_[1], y2, f);
//			y0 = term_f< 1>(w_[0], y1, f);
//			y2 = term_f<-1>(o, a_[0], y0) - shape(y1);

//		Newton's method:
//			auto const num = term_f<-1>(o, a_[0], y0) - y2 - shape<0>(y1);
//			auto const nom = term_f< 1>(1, f, term_f< 1>(shape<1>(y1), f, a_[0]);
//			y2 += num/nom;

		}
		/***/
		/*/
		template <int N_ord=0>
		XTAL_DEF_(inline)
		XTAL_LET shape(auto &y_)
		XTAL_0EX
		{
//			U_data &y0 = cache[0], &a_[0] = cache[4];
//			U_data &y1 = cache[1], &a_[1] = cache[5];
//			U_data &y2 = cache[2], &a2 = cache[6];
//			U_data &w_[0] = cache[3], &w_[1] = cache[7];
//
//			XTAL_IF0
//			XTAL_0IF (N_ord == 0) {
//				return y1*(a_[1] - 2) + antisaturator< 1, 1>(y1);
//			}
//			XTAL_0IF (N_ord == 1) {
//				return    (a_[1] - 2) + antisaturator<-1, 1>(y1);
//			}
		}
		template <int N_ord=0, int N_two=0>
		XTAL_DEF_(inline)
		XTAL_LET antisaturator(auto const &y1)
		XTAL_0EX -> void
		{
			return y1;

//			U_alpha const n_drive = S_::template head<0>();
//			U_alpha const n_curve = S_::template head<1>();
//
//			XTAL_IF0
//			XTAL_0IF (N_ord == 0) {
//				return F (y1, n_drive, n_curve)*_op::diplo_f(N_two - 1);
//			}
//			XTAL_0IF (N_ord == 1) {
//				return F'(y1, n_drive, n_curve)*_op::diplo_f(N_two - 1);
//			}
		}
		/***/

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
