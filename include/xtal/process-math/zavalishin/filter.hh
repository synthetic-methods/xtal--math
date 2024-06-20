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
:	bond::compose<filter<As>..., bond::tag<filter>>
{
};
template <typename A>
struct filter<A>
:	bond::compose<A>
{
};
template <class U_data> requires is_q<devolved_t<U_data>, U_data>
struct filter<U_data[2]>
{
	using _op = bond::operate<U_data>;
	using Z_sigma = typename _op::sigma_type;
	using Z_alpha = typename _op::alpha_type;
	static constexpr Z_alpha Z_1{1};

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

		using Z_source = _std::array<U_data, 8>;
		using Z_target = algebra::scalar_t<U_data(&)[3]>;

		Z_source cache{0, 0, 0, 0, 1, 1, 1, 0};
		Z_target y_{point_f<0>(cache), point_f<3>(cache)};
		Z_target a_{point_f<4>(cache), point_f<7>(cache)};

	public:
		using S_::S_;

		template <int N_ord=0, int N_two=0>
		XTAL_DEF_(inline)
		XTAL_LET antisaturator(U_data y1)
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
		template <int N_ord=0>
		XTAL_DEF_(inline)
		XTAL_LET shape(U_data y1)
		XTAL_0EX
		{
//			U_data &y0 = cache[0], &a0 = cache[4];
//			U_data &y1 = cache[1], &a1 = cache[5];
//			U_data &y2 = cache[2], &a2 = cache[6];
//			U_data &v0 = cache[3], &v1 = cache[7];
//
//			XTAL_IF0
//			XTAL_0IF (N_ord == 0) {
//				return y1*(a1 - 2) + antisaturator< 1, 1>(y1);
//			}
//			XTAL_0IF (N_ord == 1) {
//				return    (a1 - 2) + antisaturator<-1, 1>(y1);
//			}
		}
		template <auto ...Is>
		XTAL_DEF_(inline)
		XTAL_LET solve(U_data u, Z_alpha t)
		XTAL_0EX -> void
		{
//		With:
//			f(y[1]) = y1*(a1 - 2) + 2*F[y1, u, v]
//		
//		Fixed-point:
//			y1 = term_f< 1>(v1, y2, t);
//			y0 = term_f< 1>(v0, y1, t);
//			y2 = term_f<-1>(u, a0, y0) - shape(y1);

//		Newton's method:
//			auto const num = term_f<-1>(u, a0, y0) - y2 - shape<0>(y1);
//			auto const nom = term_f< 1>(1, t, term_f< 1>(shape<1>(y1), t, a0);
//			y2 += num/nom;

		}
		template <auto ...Is>
		XTAL_DEF_(return)
		XTAL_RET functor(U_data u, Z_alpha t)
		XTAL_0EX
		{
			using namespace horner;

			U_data &y0 = cache[0], &a0 = cache[4];
			U_data &y1 = cache[1], &a1 = cache[5];
			U_data &y2 = cache[2], &a2 = cache[6];
			U_data &v0 = cache[3], &v1 = cache[7];

			Z_alpha const t_a0 =        a0          ;// a0
			Z_alpha const t_a1 = term_f(a1, t, t_a0);// a1 + t %
			Z_alpha const t_a2 = term_f(a2, t, t_a1);// a2 + t %
			y2 = term_f<-1>(term_f<-1>((u), v0, t_a0), v1, t_a1)/t_a2;
			
			solve<Is...>(u, t);// TODO: Save current `y1`!
			
			y1 = term_f(v1, t, y2); v1 = term_f(y1, t, y2);
			y0 = term_f(v0, t, y1); v0 = term_f(y0, t, y1);

			return static_cast<Z_target const &>(y_);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_RET functor(U_data u, Z_alpha t, Z_alpha a)
		XTAL_0EX
		{
			a_[1] = a;
			return functor<Is...>(XTAL_REF_(u), t, a);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
