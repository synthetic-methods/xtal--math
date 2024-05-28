#pragma once
#include "./any.hh"
#include "../square.hh"
#include "../wrap.hh"




XTAL_ENV_(push)
namespace xtal::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` by `(-1)^(2 #) &`; spiritually equivalent to `1^# &`. \

///\param M_ism \f$\in {1, 2}\f$ specifies the underlying morphism, \
generating either circular or hyperbolic `{cosine, sine}` pairs. \

template <int M_ism=1, typename ...As> struct unity {static_assert(M_ism);};
template <int M_ism=1, typename ...As> using  unity_t = process::confined_t<unity<M_ism, As...>>;

template <int M_ism, bond::compose_q ...As> requires some_q<As...>
struct unity<M_ism, As...>: process::chain<unity<M_ism>, As...> {};


namespace _detail
{///////////////////////////////////////////////////////////////////////////////
///\
Serves as the mathematical definition of the approximant, \
which is argument-restricted by the main definition. \

template <int M_ism=1, int N_car=0> struct semiunity {static_assert(M_ism);};
template <int M_ism               > struct semiunity<M_ism,-0>: bond::compose<discard<1>, semiunity<M_ism,-1>> {};
template <int M_ism               > struct semiunity<M_ism,-1>: bond::compose<discard<2>, semiunity<M_ism,-2>> {};
template <int M_ism=1, int N_car=0> using  semiunity_t = process::confined_t<semiunity<M_ism, N_car>>;
template <int M_ism               >
struct semiunity<M_ism,-2>
{
	XTAL_LET_(int) I_sgn = sign_n<M_ism&1^1, -1>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_FN2 function(simplex_field_q auto const &w)
		XTAL_0EX
		{
			int constexpr I_lim = N_lim&0x7;

			using X = XTAL_TYP_(w); using op = bond::operate<X>;

			using alpha_t = typename op::alpha_t;
			using sigma_t = typename op::sigma_t;
			using delta_t = typename op::delta_t;

			XTAL_IF0
			XTAL_0IF_(I_lim == 0x0) {// 0:1 D[...0]@0 && D[...0]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;
				auto x = horner::polynomial_f<I_sgn>(w, x0);
				auto y = horner::polynomial_f<I_sgn>(w, y0);
				return complexion_f(x, y);
			}
			XTAL_0IF_(I_lim == 0x1) {// 2:1 D[...1]@0 && D[...0]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr x1 = 3.433629385640827046149426466881988463L;
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1);
				auto y = horner::polynomial_f<I_sgn>(w, y0);
				return complexion_f(x, y);
			}
			XTAL_0IF_(I_lim == 0x2) {// 2:3 D[...1]@0 && D[...1]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr x1 = 3.968985697854261420737444775816737932L;
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_t constexpr y1 = 2.141425248853737498352073235738997875L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1);
				return complexion_f(x, y);
			}
			XTAL_0IF_(I_lim == 0x3) {// 4:3 D[...2]@0 && D[...2]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr x1 = 4.237909660449236428693854535224051104L;
				alpha_t constexpr x2 = 0.955351046282004540229507364731116905L;
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_t constexpr y1 = 2.978283337663136395120335432185471337L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1);
				return complexion_f(x, y);
			}
			XTAL_0IF_(I_lim == 0x4) {// 4:5 D[...2]@0 && D[...3]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr x1 = 4.390338960642080575102860479411085957L;
				alpha_t constexpr x2 = 1.561426046129173147274250250239938996L;
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_t constexpr y1 = 3.457231542101901520792301595993191760L;
				alpha_t constexpr y2 = 0.331996058066891068754049734988584509L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2);
				return complexion_f(x, y);
			}
			XTAL_0IF_(I_lim == 0x5) {// 6:5 D[...3]@0 && D[...3]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr x1 = 4.487894813803155707106642350629568074L;
				alpha_t constexpr x2 = 1.975104599215236121543726210642557678L;
				alpha_t constexpr x3 = 0.094096528134246695768752645930262507L;
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_t constexpr y1 = 3.763711817027786915281807613442402050L;
				alpha_t constexpr y2 = 0.623295934882939155479821936321002832L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2, x3);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2);
				return complexion_f(x, y);
			}
			XTAL_0IF_(I_lim == 0x6) {// 6:7 D[...3]@0 && D[...4]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr x1 = 4.556288946308768410678019439942678932L;
				alpha_t constexpr x2 = 2.275378023833944536636939235400154936L;
				alpha_t constexpr x3 = 0.200888908030920297664420149309522799L;
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_t constexpr y1 = 3.978578321256066662987441200793155884L;
				alpha_t constexpr y2 = 0.859750287076040352418274981236440341L;
				alpha_t constexpr y3 = 0.022706582386768037080226242510354953L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2, x3);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2, y3);
				return complexion_f(x, y);
			}
			XTAL_0IF_(I_lim == 0x7) {// 8:7 D[...4]@0 && D[...5]@¼
				alpha_t constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_t constexpr x1 = 4.606547402370024772724248519073301565L;
				alpha_t constexpr x2 = 2.500953736670776238483254500604284293L;
				alpha_t constexpr x3 = 0.300243870518009486091801923099845388L;
				alpha_t constexpr x4 = 0.004749448374844164315053279863807317L;
				alpha_t constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_t constexpr y1 = 4.136469917598875065501883380414805209L;
				alpha_t constexpr y2 = 1.048974757862995076084471175276645114L;
				alpha_t constexpr y3 = 0.054095324243024101847331195900039153L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2, x3, x4);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2, y3);
				return complexion_f(x, y);
			}
		}

	};
};


}///////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct unity<M_ism>
{
//	XTAL_LET_(int) I_sgn = sign_n<M_ism&1^1, -1>;

	using subkind = process::chain<square<M_ism, 0>, _detail::semiunity<M_ism,-0>>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_FN2 function(complex_field_q auto const &t)
		XTAL_0EX
		{
			return function<N_lim>(t.real(), t.imag());
		}
		template <int N_lim=-1>
		XTAL_FN2 function(auto &&t_1, simplex_field_q auto &&t_i)
		XTAL_0EX
		{
			using _std::exp;

			using T_i = XTAL_TYP_(t_i); using op = bond::operate<T_i>;
			return function<N_lim>(XTAL_REF_(t_1))*exp(XTAL_REF_(t_i)*op::patio_f(-2));
		}

		template <int N_lim=-1>
		XTAL_FN2 function(simplex_field_q auto o)
		XTAL_0EX
		{
			using Op = XTAL_TYP_(o); using op = bond::operate<Op>;

			if constexpr (N_lim < 0) {
				using namespace _std;

				auto w = o*op::patio_2;
				XTAL_IF0
				XTAL_0IF_(1 == M_ism) {return complexion_f(cos (w), sin (w));}
				XTAL_0IF_(2 == M_ism) {return complexion_f(cosh(w), sinh(w));}
			}
			else {
				XTAL_VAL f_assign = [] XTAL_1FN_(op::assign_f);
				XTAL_USE F_assign = decltype(f_assign);
				
				auto w = wrap_f(o);
				auto m = wrap_f(w*op::diplo_1)*op::haplo_1;
				return S_::template function<N_lim>(m)*hoist_f(f_assign, m != w);
			}
		}
		template <int N_lim=-1>
		XTAL_FN2 function(algebra::differential::circular_q auto w)
		XTAL_0EX
		{
			using Op = XTAL_TYP_(w); using op = bond::operate<Op>;
			XTAL_USE T_alpha = typename op::alpha_t;
			XTAL_USE F_alpha = decltype([] XTAL_1FN_(_std::bit_cast<T_alpha>));

			if constexpr (N_lim < 0) {
				return function<N_lim>(w(0));
			}
			else {
				auto const bit_sign = w[0] & op::sign.mask;
				w[0] &= op::positive.mask >> 1;
				w[0] |=          bit_sign >> 1;
				return S_::template function<N_lim>(w(0))*hoist_f<F_alpha>(bit_sign|op::unit.mask);
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
