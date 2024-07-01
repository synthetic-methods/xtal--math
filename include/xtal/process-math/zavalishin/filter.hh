#pragma once
#include "./any.hh"

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
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	//	NOTE: Expected maximum is 64/8: 6 doubles not including coefficients...
		XTAL_SET N_cache = bond::operating::alignment{}*sizeof(size_type);
		alignas (N_cache) _std::byte m_cache[N_cache];

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return)
		XTAL_LET method(auto &&o, auto &&f)
		XTAL_0EX
		{
			using namespace horner;
			using V  = XTAL_ALL_(o); using _op = bond::operate<V>;
			using A  = typename bond::operate<V>::alpha_type;
			
			using Y_ = algebra::scalar_t<V[3]>;
			using V_ = algebra::scalar_t<V[2]>;
			using A_ = algebra::scalar_t<A[2]>;
			
			static_assert(aligned_n<V_, A_> <= N_cache);
			static_assert(_std::is_trivially_destructible_v<V>);
			
			size_type i{};
			V_ &v_ = reinterpret_cast<V_ &>(m_cache[maligned_f<V_>(i)]);
			A_ &a_ = reinterpret_cast<V_ &>(m_cache[maligned_f<A_>(i)]);

			A const &a0 = a_[0];
			A const &a1 = a_[1];
			A const  a2 =    1 ;

			A const w0 =        a0        ;// a_[0]
			A const w1 = term_f(a1, w0, f);// a_[1] + f %
			A const w2 = term_f(a2, w1, f);// a_[2] + f %
			Y_ y_{V{}, V{}, term_f<-1>(term_f<-1>(o, v_[0], w0), v_[1], w1)/w2};

		//	solve<Is...>(y_, o, f);// TODO: Save current `y_[1]`!
			
			y_[1] = term_f(v_[1], y_[2], f); v_[1] = term_f(y_[1], y_[2], f);
			y_[0] = term_f(v_[0], y_[1], f); v_[0] = term_f(y_[0], y_[1], f);

			return y_;
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&o, auto &&f, auto &&_R)
		XTAL_0EX
		{
			using V  = XTAL_ALL_(o); using _op = bond::operate<V>;
			using A  = typename _op::alpha_type;

			using V_ = algebra::scalar_t<V[2]>;
			using A_ = algebra::scalar_t<A[2]>;

			size_type i{};
			V_   &v_ = reinterpret_cast<V_ &>(m_cache[maligned_f<V_>(i)]);
			A_   &a_ = reinterpret_cast<V_ &>(m_cache[maligned_f<A_>(i)]);
			new (&a_) A_{_op::alpha_1, _op::alpha_2*_R};
			return method<Is...>(XTAL_REF_(o), f);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&o, auto &&f, complex_number_q auto &&_S)
		XTAL_0EX
		{
			//\todo\
			Use `differential::circular` to provide the base frequency, \
			combining relative-frequency and resonance as the complex `s`, \
			where frequency is given by `Abs[s]`, and `Q = Sqrt[# - ¼]@?`. \

			using V  = XTAL_ALL_(o); using _op = bond::operate<V>;
			using A  = typename _op::alpha_type;

			using _std::exp;// Can approximate?
			auto constexpr N_lnH = (A) -0.6931471805599453094172321214581765681e+0L;
			
			return method(XTAL_REF_(o), XTAL_REF_(f), exp(N_lnH*XTAL_REF_(_S).imag()));
		}
		template <auto ...Is>
		XTAL_DEF_(inline)
		XTAL_LET solve(auto &y_, auto const &o, auto const &f)
		XTAL_0EX -> void
		{
//		With:
//			f(y[1]) = y1*(a_[1] - 2) + 2*F[y1, o, v]
//		
//		Fixed-point:
//			y1 = term_f< 1>(v_[1], y2, f);
//			y0 = term_f< 1>(v_[0], y1, f);
//			y2 = term_f<-1>(o, a_[0], y0) - shape(y1);

//		Newton's method:
//			auto const num = term_f<-1>(o, a_[0], y0) - y2 - shape<0>(y1);
//			auto const nom = term_f< 1>(1, f, term_f< 1>(shape<1>(y1), f, a_[0]);
//			y2 += num/nom;

		}
		template <int N_ord=0>
		XTAL_DEF_(inline)
		XTAL_LET shape(auto &y_)
		XTAL_0EX
		{
//			U_data &y0 = cache[0], &a_[0] = cache[4];
//			U_data &y1 = cache[1], &a_[1] = cache[5];
//			U_data &y2 = cache[2], &a2 = cache[6];
//			U_data &v_[0] = cache[3], &v_[1] = cache[7];
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

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
