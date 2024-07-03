#pragma once
#include "./any.hh"
#include "../square.hh"
#include "../wrap.hh"
#include "../../occur-math/all.hh"



XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` by `(-1)^(2 #) &`; spiritually equivalent to `1^# &`. \

///\param M_ism \f$\in {1, 2}\f$ specifies the underlying morphism, \
generating either circular or hyperbolic `{cosine, sine}` pairs. \

template <int M_ism=1, typename ...As> XTAL_TYP unity {static_assert(M_ism);};
template <int M_ism=1, typename ...As> XTAL_USE unity_t = process::confined_t<unity<M_ism, As...>>;

using V_unity_limit = occur::math::limit_t<2>;


namespace _detail
{///////////////////////////////////////////////////////////////////////////////
///\
Serves as the mathematical definition of the approximant, \
which is argument-restricted by the main definition. \

template <int M_ism=1, int N_car=0> XTAL_TYP subunity {static_assert(M_ism);};
template <int M_ism               > XTAL_TYP subunity<M_ism,-0>: bond::compose<discarding<1, +1>, subunity<M_ism,-1>> {};
template <int M_ism               > XTAL_TYP subunity<M_ism,-1>: bond::compose<discarding<1, +2>, subunity<M_ism,-2>> {};
template <int M_ism=1, int N_car=0> XTAL_USE subunity_t = process::confined_t<subunity<M_ism, N_car>>;
template <int M_ism               >
struct subunity<M_ism,-2>
{
	static constexpr int I_sgn = sign_n<(M_ism&1)^1, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,static)
		XTAL_LET function(simplex_field_q auto const &w)
		XTAL_0EX -> decltype(auto)
		{
			int constexpr I_lim = N_lim&0x7;

			using X = XTAL_ALL_(w); using _op = bond::operate<X>;

			using alpha_type = typename _op::alpha_type;
			using sigma_type = typename _op::sigma_type;
			using delta_type = typename _op::delta_type;

			XTAL_IF0
			XTAL_0IF (I_lim == 0x0) {// 0:1 D[...0]@0 && D[...0]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;
				auto x = horner::polynomial_f<I_sgn>(w, x0);
				auto y = horner::polynomial_f<I_sgn>(w, y0);
				return complexion_f(x, y);
			}
			XTAL_0IF (I_lim == 0x1) {// 2:1 D[...1]@0 && D[...0]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr x1 = 3.433629385640827046149426466881988463L;
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1);
				auto y = horner::polynomial_f<I_sgn>(w, y0);
				return complexion_f(x, y);
			}
			XTAL_0IF (I_lim == 0x2) {// 2:3 D[...1]@0 && D[...1]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr x1 = 3.968985697854261420737444775816737932L;
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_type constexpr y1 = 2.141425248853737498352073235738997875L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1);
				return complexion_f(x, y);
			}
			XTAL_0IF (I_lim == 0x3) {// 4:3 D[...2]@0 && D[...2]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr x1 = 4.237909660449236428693854535224051104L;
				alpha_type constexpr x2 = 0.955351046282004540229507364731116905L;
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_type constexpr y1 = 2.978283337663136395120335432185471337L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1);
				return complexion_f(x, y);
			}
			XTAL_0IF (I_lim == 0x4) {// 4:5 D[...2]@0 && D[...3]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr x1 = 4.390338960642080575102860479411085957L;
				alpha_type constexpr x2 = 1.561426046129173147274250250239938996L;
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_type constexpr y1 = 3.457231542101901520792301595993191760L;
				alpha_type constexpr y2 = 0.331996058066891068754049734988584509L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2);
				return complexion_f(x, y);
			}
			XTAL_0IF (I_lim == 0x5) {// 6:5 D[...3]@0 && D[...3]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr x1 = 4.487894813803155707106642350629568074L;
				alpha_type constexpr x2 = 1.975104599215236121543726210642557678L;
				alpha_type constexpr x3 = 0.094096528134246695768752645930262507L;
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_type constexpr y1 = 3.763711817027786915281807613442402050L;
				alpha_type constexpr y2 = 0.623295934882939155479821936321002832L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2, x3);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2);
				return complexion_f(x, y);
			}
			XTAL_0IF (I_lim == 0x6) {// 6:7 D[...3]@0 && D[...4]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr x1 = 4.556288946308768410678019439942678932L;
				alpha_type constexpr x2 = 2.275378023833944536636939235400154936L;
				alpha_type constexpr x3 = 0.200888908030920297664420149309522799L;
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_type constexpr y1 = 3.978578321256066662987441200793155884L;
				alpha_type constexpr y2 = 0.859750287076040352418274981236440341L;
				alpha_type constexpr y3 = 0.022706582386768037080226242510354953L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2, x3);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2, y3);
				return complexion_f(x, y);
			}
			XTAL_0IF (I_lim == 0x7) {// 8:7 D[...4]@0 && D[...5]@¼
				alpha_type constexpr x0 = 1.000000000000000000000000000000000000L;// 1
				alpha_type constexpr x1 = 4.606547402370024772724248519073301565L;
				alpha_type constexpr x2 = 2.500953736670776238483254500604284293L;
				alpha_type constexpr x3 = 0.300243870518009486091801923099845388L;
				alpha_type constexpr x4 = 0.004749448374844164315053279863807317L;
				alpha_type constexpr y0 = 3.141592653589793238462643383279502884L;// pi
				alpha_type constexpr y1 = 4.136469917598875065501883380414805209L;
				alpha_type constexpr y2 = 1.048974757862995076084471175276645114L;
				alpha_type constexpr y3 = 0.054095324243024101847331195900039153L;
				auto x = horner::polynomial_f<I_sgn>(w, x0, x1, x2, x3, x4);
				auto y = horner::polynomial_f<I_sgn>(w, y0, y1, y2, y3);
				return complexion_f(x, y);
			}
		}

	};
};


}///////////////////////////////////////////////////////////////////////////////

template <int M_ism, bond::compose_q ...As> requires some_q<As...>
struct unity<M_ism, As...>
:	process::lift<unity<M_ism>, bond::compose<As...>>
{
};
template <int M_ism> requires (0 < M_ism)
struct unity<M_ism>
{
	using superprocess = process::lift_t<square<M_ism, 0>, _detail::subunity<M_ism,-0>>;

	using subkind = bond::compose<void
	,	V_unity_limit::dispatch<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=-1, class U>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(_std::initializer_list<U> o)
		XTAL_0EX -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			_std::complex<U> w; auto &m = involved_f(w);
			_std::copy_n(point_f(o), 2, m);
			return function<N_lim>(w);
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(complex_field_q auto const &t)
		XTAL_0EX -> decltype(auto)
		{
			return function<N_lim>(t.real(), t.imag());
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&t_re, simplex_field_q auto &&t_im)
		XTAL_0EX -> decltype(auto)
		{
			using _std::exp;

			using _op = bond::operate<decltype(t_re), decltype(t_im)>;
			return function<N_lim>(XTAL_REF_(t_re))*exp(XTAL_REF_(t_im)*_op::patio_f(-2));
		}

		template <int N_lim=-1>
		XTAL_DEF_(return,static)
		XTAL_LET function(simplex_field_q auto o)
		XTAL_0EX -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			using Op_alpha = typename _op::alpha_type;

			if constexpr (N_lim < 0) {
				using namespace _std;

				auto w = o*_op::patio_2;
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return complexion_f(cos (w), sin (w));}
				XTAL_0IF (2 == M_ism) {return complexion_f(cosh(w), sinh(w));}
			}
			else {
				auto constexpr assigned_f = [] (int i) XTAL_0FN -> Op_alpha {return (i << 1) - 1;};
				auto w = wrap_f(o);
				auto m = wrap_f(w*_op::diplo_1)*_op::haplo_1;
				return superprocess::template function<N_lim>(m)*operative_f(assigned_f, m == w);
			}
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,static)
		XTAL_LET function(algebra::d_::circular_q auto d)
		XTAL_0EX -> decltype(auto)
		{
			using _op = bond::operate<decltype(d)>;
			XTAL_USE Op_alpha = typename _op::alpha_type;
			XTAL_USE Fn_alpha = decltype([] XTAL_1FN_(_xtd::bit_cast<Op_alpha>));

			if constexpr (N_lim < 0) {
				return function<N_lim>(d(0));
			}
			else {
				auto &d0 = d[0];
				auto  s0 = d0 << 0;
				auto  s1 = d0 << 1;
				auto  sn = s0^s1; sn &= _op::sign.mask;
				d0 ^= sn;
				return superprocess::template function<N_lim>(d(0))*inoperative_f<Fn_alpha>(_op::unit.mask|sn);
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
