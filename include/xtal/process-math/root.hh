#pragma once
#include "./any.hh"

#include "./cut.hh"
#include "./magnum.hh"
#include "./square.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_exp=1, int M_cut=0>
struct   root;

template <int M_exp=1, int M_cut=0>
using    root_t = process::confined_t<root<M_exp, M_cut>>;

template <int M_exp=1, int M_cut=0, auto ...Ns>
XTAL_DEF_(short)
XTAL_LET root_f(auto &&w)
noexcept -> decltype(auto)
{
	XTAL_IF0
	XTAL_0IF (M_exp == 1) {return                                                XTAL_REF_(w) ;}
	XTAL_0IF (M_exp != 1) {return root_t<M_exp, M_cut>::template function<Ns...>(XTAL_REF_(w));}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_exp, int M_cut> requires (0 <  M_cut)
struct root<M_exp, M_cut>
{
	using superkind = root<M_exp>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&w)
		noexcept -> auto
		{
			using _op = bond::operate<XTAL_ALL_(w)>;
			XTAL_LET N_exp = magnum_f(M_exp*M_cut);
			XTAL_LET N_min = _op::minilon_f(N_exp);
			return S_::template function<Ns...>(cut_f<XTAL_VAL_(N_min)>(XTAL_REF_(w)));
		}

	};
};
template <int M_exp, int M_cut> requires (M_cut <= 0)
struct root<M_exp, M_cut>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0b11> requires un_q<M_exp, -3>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&w)
		noexcept -> auto
		{
			using W_op = bond::operate<XTAL_ALL_(w)>;
			XTAL_IF0
			XTAL_0IF (integral_q<decltype(w)>) {
				return function<N_lim>(W_op::alpha_f(w));
			}
			XTAL_0IF (M_exp ==  1) {return               XTAL_REF_(w)  ;}
			XTAL_0IF (M_exp ==  2) {return          sqrt(XTAL_REF_(w)) ;}
			XTAL_0IF (M_exp ==  3) {return          cbrt(XTAL_REF_(w)) ;}
			XTAL_0IF (M_exp ==  4) {return     sqrt(sqrt(XTAL_REF_(w)));}
			XTAL_0IF (M_exp == -1) {return one/          XTAL_REF_(w)  ;}
			XTAL_0IF (M_exp == -2) {return one/     sqrt(XTAL_REF_(w)) ;}
			XTAL_0IF (M_exp == -3) {return one/     cbrt(XTAL_REF_(w)) ;}
			/**/
			XTAL_0IF (M_exp == -4) {return one/sqrt(sqrt(XTAL_REF_(w)));}
			/*/
			XTAL_0IF (M_exp == -4) {
				auto constexpr c0 =   W_op::ratio_f( 5, 4);
				auto const     c1 = w*W_op::ratio_f(-1, 4);
				auto y = XTAL_REF_(w);
				y *= W_op::template root_f<-2>(y);
				y  = _op::template root_f<-2>(y);
				y *= _xtd::accumulator(c0, c1, square_f(square_f(y)));
				y *= _xtd::accumulator(c0, c1, square_f(square_f(y)));
				return y;
			}
			/***/
			XTAL_0IF (0 == M_exp%2) {
				return root_f<M_exp/2>(root_f<2>(XTAL_REF_(w)));
			}
			XTAL_0IF (0 == M_exp%3) {
				return root_f<M_exp/3>(root_f<3>(XTAL_REF_(w)));
			}
			XTAL_0IF_(else) {
				return pow(XTAL_REF_(w), W_op::alpha_1/M_exp);
			}
		}
		template <int N_lim=0b11> requires in_q<M_exp, -3>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&w)
		noexcept -> auto
		{
			using W    = XTAL_ALL_(w);
			using W_op = bond::operate<W>;

			XTAL_IF0
			XTAL_0IF (integral_q<decltype(w)>) {
				return function<N_lim>(W_op::alpha_f(w));
			}
			XTAL_0IF (N_lim < 0 or not real_number_q<decltype(w)>) {
				return one/cbrt(XTAL_REF_(w));
			}
			XTAL_0IF_(else) {
				using X_delta = typename W_op::delta_type;
				using X_sigma = typename W_op::sigma_type;
				using X_alpha = typename W_op::alpha_type;

				X_sigma constexpr N     =   W_op::sign.mask >> 11;    // {64,32} -> {52,20}
				X_sigma constexpr K_exp = 3*W_op::full.depth/4 - 1; // {64,32} -> {47,23}
				X_sigma constexpr K_num = 37;
				X_sigma constexpr K_nom = 30;
				X_sigma constexpr K     = (K_num << K_exp)/K_nom;// 4/3 - 1/10 // Together
			//	X_sigma constexpr K     = 0x9DDDDD;// 32-bit
			//	X_sigma constexpr K     = 0x9DDDDDDDDDDD;// 64-bit

			//	X_sigma constexpr magic = 0x3FF*(1/3)*(3 - 1)*(5 - 3)*N - K;
				X_sigma constexpr magic = 0x554*N - K;

				auto m = _xtd::bit_cast<X_sigma>(w);
				auto v = W_op::sign.mask&m;
				m ^= v;
				m /= 3;
				m  = magic - m;
				m ^= v;
				auto y = _xtd::bit_cast<X_alpha>(m);
				
				auto constexpr a = W_op::ratio_f( 4, 3);
				auto const     x = W_op::ratio_f(-1, 3)*XTAL_REF_(w);
				auto constexpr I = below_m<(1<<2), (unsigned) N_lim>;
				#pragma unroll
				for (unsigned i{}; i < I; ++i) {
					y *= _xtd::accumulator(a, x, y*y*y);
				}
				return y;
			}
		}
		/*/
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_number_q auto &&w)
		noexcept -> auto
		{
			using _op = bond::operate<XTAL_ALL_(w)>;
			auto constexpr N_sqrt_half = (typename _op::alpha_type) 0.7071067811865475244008443621048490393L;

			XTAL_IF0
			XTAL_0IF (M_exp ==  1) {
				return XTAL_REF_(w);
			}
			XTAL_0IF (M_exp ==  2) {
				auto const x   = w.real();
				auto const y   = w.imag();
				auto const n   = function<Ns...>(x*x + y*y);
				auto const lhs = function<Ns...>(n + x);
				auto const rhs = function<Ns...>(n - x);
				return N_sqrt_half*complexion_f(lhs, rhs*_op::assigned_f(y));
			}
			XTAL_0IF (M_exp == -1) {
				return one/XTAL_REF_(w);
			}
			XTAL_0IF (M_exp == -2) {
				using dis = root_t<2>;
				auto       x   =  w.real();
				auto       y   =  w.imag();
				auto const a    =  x*x + y*y;
				auto const u_dn =  function<Ns...>(a);
				auto const u_up =  u_dn*a;
				auto const v_re =  u_dn*dis::template function<Ns...>(u_up + x);
				auto const v_im = -u_dn*dis::template function<Ns...>(u_up - x);

				return N_sqrt_half*complexion_f(v_re, v_im*_op::assigned_f(y));
			}
		}
		/***/

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
