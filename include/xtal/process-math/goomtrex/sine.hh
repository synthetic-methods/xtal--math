#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::goomtrex
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/***
Evaluates the hyperbolic-sine approximation:
```
Zoom[Inv@Sqrt[(n - 1) (n + 1)]][#*1&, # + Sqrt@Term[1, #, #]&, #^n^+1&, Half[# - #^-1] &, #/n&]
Zoom[Inv@Sqrt[(n - 1) (n + 1)]][#*n&, # + Sqrt@Term[1, #, #]&, #^n^-1&, Half[# - #^-1] &, #/1&]
```
***/
template <auto  ...Ms>	struct   sine;
template <auto  ...Ms>	using    sine_t = process::confined_t<sine<Ms...>>;

////////////////////////////////////////////////////////////////////////////////

template <auto  ...Ms>
struct   sine
{
	template <class S>
	using subtype = bond::compose_s<S>;

};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, int M_car> requires in_q<M_ism, 1, 2> and in_q<M_car, -0>
struct   sine<M_ism, M_car>
:	bond::compose<discarded<1>, sine<M_ism, -1>>
{
};
template <int M_ism, int M_car> requires in_q<M_ism, 1, 2> and in_q<M_car, -1>
struct   sine<M_ism, M_car>
{
	using supertype = bond::compose<discarded<2, M_ism>, sine<M_ism, -2>>;

	template <class S>
	class subtype : public bond::compose_s<S, supertype>
	{
		using S_ = bond::compose_s<S, supertype>;

	public:
		using S_::S_;

		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto u)
		noexcept -> auto
		{
			using  U_op = bond::operate<decltype(u)>;
			XTAL_LET      a0 = U_op::alpha_1;
			XTAL_LET haplo_0 = U_op::haplo_0;
			XTAL_LET haplo_1 = U_op::haplo_1;

			XTAL_LET  N     = magnum_n<N_ord>;
			XTAL_LET  N_par = N_ord&1;
			XTAL_IF0
			XTAL_0IF (N_ord <  0) {
				return sine_t<M_ism, M_car + 1>::template function<N_ord>(u)*root_f<-1, 1>(u);
			}
			XTAL_0IF (N_par == 0) {
				auto constexpr zoom  = U_op::alpha_f((N_ord - 1)*(N_ord + 1));
				auto constexpr zoom_ = U_op::template roots_f<2>(zoom);
				auto const [zoom_dn, zoom_up] = zoom_;
				if (N_ord < 0) {u *= U_op::ratio_f(N, 1);}
				u *=   zoom_up;
				u +=   root_f< 2>(term_f(U_op::alpha_1, u, u));
				u  = nomial_f< N>(u);
				u -=   root_f<-1>(u);
				u *=   haplo_1;
				u *=   zoom_dn;
				if (0 < N_ord) {u *= U_op::ratio_f(1, N);}
				return u;
			}
			XTAL_0IF (N_ord == 9) {
			//	Module[{u=#, w=# #}, Module[{t0=Term[1, 1/60, w]}, Module[{t1=u*t0*Sqrt[3/20]}, Term[1, t1, t1]*t0]]]&
				auto const w = square_f(u);
				XTAL_LET  a1 = U_op::alpha_f(0.1666666666666666666666666666666667e-1L);// 1/60
				XTAL_LET  b1 = U_op::alpha_f(0.3872983346207416885179265399782400e-0L);// Sqrt[3/20]

				auto const t0 = term_f(a0, a1, w), tX = t0*b1*u;
				return t0*term_f(a0, tX, tX);
			}
//			XTAL_0IF (N_ord == 27) {
//				auto constexpr a = U{0.2336664289109584522132436978521609e2L};// Sqrt[1/546]
//				auto const o = a * u;
//				auto const t = term_f(1, o, o);
//				auto const s = square_f(t*o*U_co(9, 1));
//				auto const r = term_f(1, s, U_co(1, 9));
//				auto const y = term_f(1, s, r*r)*r*t;
//				return y;
//			}
			XTAL_0IF_(else) {
				return S_::template function<N_ord>(XTAL_MOV_(u));
			}
		}
		template <int N_ord=0> requires (N_ord <  0)
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> auto
		{
			return sine_t<-M_ism, M_car>::template function<-N_ord>(XTAL_REF_(u));
		}
		template <int N_ord=0> requires (N_ord == 0)
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> auto
		{
			return XTAL_REF_(u);
		}

	};
};
template <int M_ism, int M_car> requires in_q<M_ism, 1, 2> and in_q<M_car, -2>
struct   sine<M_ism, M_car>
{
	using supertype = sine<M_ism, -4>;

	template <class S>
	class subtype : public bond::compose_s<S, supertype>
	{
		using S_ = bond::compose_s<S, supertype>;

	public:
		using S_::S_;

		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto w)
		noexcept -> decltype(auto)
		{
			using  W_op = bond::operate<decltype(w)>;
			XTAL_LET a0 = W_op::alpha_1;

			XTAL_LET  N_par = N_ord&1;
			XTAL_IF0
			XTAL_0IF (0 == N_par or N_ord < 0) {
				auto const u = root_f<2>(magnum_f(XTAL_MOV_(w)));
				return sine_t<M_ism, M_car + 1>::template function<N_ord>(u);
			}
			XTAL_0IF (2 == N_ord) {
				XTAL_LET a1 = W_op::ratio_f(1, 3);
				return root_f<2>(term_f(a0, a1, w));
			}
			XTAL_0IF (3 == N_ord) {
				XTAL_LET a1 = W_op::ratio_f(1, 6);
				return          (term_f(a0, a1, w));
			}
			XTAL_0IF (4 == N_ord) {
				XTAL_LET a1 = W_op::ratio_f(1, 15);
				XTAL_LET a2 = W_op::ratio_f(2, 15);
				return root_f<2>(term_f(a0, a1, w))*term_f(a0, w, a2);
			}
		//	XTAL_0IF (5 == N_ord) {
		//		XTAL_LET a0 = W_op::alpha_f(  5.21298012204885);
		//		XTAL_LET a1 = W_op::alpha_f(207.17516175287642);
		//		return term_f(a1, w, w - a0)*(w + a0);
		//	}
		//	XTAL_0IF (7 == N_ord) {
		//		XTAL_LET a0 = W_op::alpha_f( 9.036244755390396);
		//		XTAL_LET a1 = W_op::alpha_f(29.340502414951544);
		//		XTAL_LET a2 = W_op::alpha_f(45.623252829658050);
		//		return (w + a0)*(w + a1)*(w + a2);
		//	}
		//	XTAL_0IF (9 == N_ord) {
		//		XTAL_LET a0 = W_op::alpha_f( 9.358222275240875);
		//		XTAL_LET a1 = W_op::alpha_f(33.054072893322780);
		//		XTAL_LET a2 = W_op::alpha_f(60.000000000000000);
		//		XTAL_LET a3 = W_op::alpha_f(77.587704831436300);
		//		return (w + a0)*(w + a1)*(w + a2)*(w + a3);
		//	}
			XTAL_0IF (9 == N_ord) {
				XTAL_LET   a1 = W_op::alpha_f(0.1666666666666666666666666666666667e-1L);// 1/60
				XTAL_LET   b1 = W_op::alpha_f(0.3872983346207416885179265399782400e-0L);// Sqrt[3/20]
				auto const t0 = term_f(a0, a1, w), t1 = t0*b1;
				return t0*term_f(a0, t1, t1, w);
			}
			XTAL_0IF_(else) {
				return term_f(a0, w, S_::template function<N_ord>(w));
			}
		}
		template <int N_ord=0> requires (N_ord <  0)
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> auto
		{
			return sine_t<-M_ism, M_car>::template function<-N_ord>(XTAL_REF_(u));
		}
		template <int N_ord=0> requires (N_ord == 0)
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> auto
		{
			using U = XTAL_ALL_(u);
			return U{bond::operate<U>::alpha_1};
		}

	};
};
template <int M_ism, int M_car> requires in_q<M_ism, 1, 2> and in_q<M_car, -3>
struct   sine<M_ism, M_car>
{
	using supertype = sine<M_ism, -4>;

	template <class S>
	class subtype : public bond::compose_s<S, supertype>
	{
		using S_ = bond::compose_s<S, supertype>;

	public:
		using S_::S_;

		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto w)
		noexcept -> auto
		{
			using W_op = bond::operate<decltype(w)>;

			XTAL_LET  N_par = N_ord&1;
			XTAL_IF0
			XTAL_0IF (0 == N_par or N_ord < 0) {
				return sine_t<M_ism, M_car + 1>::template function<N_ord>(w) - W_op::alpha_1;
			}
			XTAL_0IF_(else) {
				return S_::template function<N_ord>(w)*w;
			}
		}

	};
};
template <int M_ism, int M_car> requires in_q<M_ism, 1, 2> and in_q<M_car, -4>
struct   sine<M_ism, M_car>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto w)
		noexcept -> auto
		{
			using W_op = bond::operate<decltype(w)>;

			XTAL_LET  N_par = N_ord&1;
			XTAL_IF0
			XTAL_0IF (0 == N_par or N_ord < 0) {
				return sine_t<M_ism, M_car + 1>::template function<N_ord>(w)*root_f<-1, 1>(w);
			}
			XTAL_0IF (1 == N_ord) {
				return W_op::alpha_0;
			}
			XTAL_0IF (3 == N_ord) {
				return W_op::ratio_f(1, 6);
			}
			XTAL_0IF (5 == N_ord) {
				XTAL_LET a1 = W_op::ratio_f(1,       6);
				XTAL_LET a3 = W_op::ratio_f(1,    1080);
				return term_f(a1, a3, square_f(w));
			}
			XTAL_0IF (7 == N_ord) {
				XTAL_LET a1 = W_op::ratio_f(1,       6);
				XTAL_LET a2 = W_op::ratio_f(1,     144);
				XTAL_LET a3 = W_op::ratio_f(1,   12096);
				return termial_f(w, a1, a2, a3);
			}
			XTAL_0IF (9 == N_ord) {
				XTAL_LET a1 = W_op::ratio_f(1,       6);
				XTAL_LET a2 = W_op::ratio_f(3,     400);
				XTAL_LET a3 = W_op::ratio_f(1,    8000);
				XTAL_LET a4 = W_op::ratio_f(1, 1440000);
				return termial_f(XTAL_REF_(w), a1, a2, a3, a4);
			}
			XTAL_0IF_(void)
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
/*/
template <auto  ...Ms            >
struct   sine<bond::seek_t<Ms...>{}> :	sine<Ms...> {};
/***/
template <auto     Ms, auto ...Ns>
XTAL_DEF sine_f(auto &&...oo) noexcept {return sine_t<Ms>::template function<Ns...>(XTAL_REF_(oo)...);}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
