#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Provides the argument-restricted implementations for `unity` and `tangy`. \


template <int M_ism=1, int N_car=0> struct   impunity;
template <int M_ism=1, int N_car=0> using    impunity_t = process::confined_t<impunity<M_ism>>;

template <int M_ism=1, auto ...Ns>
XTAL_DEF_(short)
XTAL_LET impunity_f(auto &&o)
noexcept -> decltype(auto)
{
	return impunity_t<M_ism>::template function<Ns...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> struct impunity<M_ism,-0>: bond::compose<discarded<1       >, impunity<M_ism,-1>> {};
template <int M_ism> struct impunity<M_ism,-1>: bond::compose<discarded<2, M_ism>, impunity<M_ism,-2>> {};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 1, 2>
struct impunity<M_ism,-2>
{
	//\
	static constexpr int I_sgn = sign_n<(M_ism&1)^1, -1>;
	static constexpr int I_sgn = 1;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(long,static)
		XTAL_LET function(simplex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using X = XTAL_ALL_(o); using _fix = bond::fixture<X>;

			using alpha_type = typename _fix::alpha_type;
			using sigma_type = typename _fix::sigma_type;
			using delta_type = typename _fix::delta_type;

			auto const w = objective_f(o);

			XTAL_IF0
			XTAL_0IF (N_lim == 0x0) {// 0:1 D[...0]@0 && D[...0]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				auto const x = termial_f(w, x0);
				auto const y = termial_f(w, y0);
				return complexion_f(x, y);
			}
			XTAL_0IF (N_lim == 0x1) {// 2:1 D[...1]@0 && D[...0]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr x1 = 3.433629385640827046149426466881988463L;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				auto const x = termial_f(w, x0, x1);
				auto const y = termial_f(w, y0);
				return complexion_f(x, y);
			}
			XTAL_0IF (N_lim == 0x2) {// 2:3 D[...1]@0 && D[...1]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr x1 = 3.968985697854261420737444775816737932L;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				alpha_type constexpr y1 = 2.141425248853737498352073235738997875L;
				auto const x = termial_f(w, x0, x1);
				auto const y = termial_f(w, y0, y1);
				return complexion_f(x, y);
			}
			XTAL_0IF (N_lim == 0x3) {// 4:3 D[...2]@0 && D[...2]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr x1 = 4.237909660449236428693854535224051104L;
				alpha_type constexpr x2 = 0.955351046282004540229507364731116905L;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				alpha_type constexpr y1 = 2.978283337663136395120335432185471337L;
				auto const x = termial_f(w, x0, x1, x2);
				auto const y = termial_f(w, y0, y1);
				return complexion_f(x, y);
			}
			XTAL_0IF (N_lim == 0x4) {// 4:5 D[...2]@0 && D[...3]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr x1 = 4.390338960642080575102860479411085957L;
				alpha_type constexpr x2 = 1.561426046129173147274250250239938996L;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				alpha_type constexpr y1 = 3.457231542101901520792301595993191760L;
				alpha_type constexpr y2 = 0.331996058066891068754049734988584509L;
				auto const x = termial_f(w, x0, x1, x2);
				auto const y = termial_f(w, y0, y1, y2);
				return complexion_f(x, y);
			}
			XTAL_0IF (N_lim == 0x5) {// 6:5 D[...3]@0 && D[...3]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr x1 = 4.487894813803155707106642350629568074L;
				alpha_type constexpr x2 = 1.975104599215236121543726210642557678L;
				alpha_type constexpr x3 = 0.094096528134246695768752645930262507L;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				alpha_type constexpr y1 = 3.763711817027786915281807613442402050L;
				alpha_type constexpr y2 = 0.623295934882939155479821936321002832L;
				auto const x = termial_f(w, x0, x1, x2, x3);
				auto const y = termial_f(w, y0, y1, y2);
				return complexion_f(x, y);
			}
			XTAL_0IF (N_lim == 0x6) {// 6:7 D[...3]@0 && D[...4]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr x1 = 4.556288946308768410678019439942678932L;
				alpha_type constexpr x2 = 2.275378023833944536636939235400154936L;
				alpha_type constexpr x3 = 0.200888908030920297664420149309522799L;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				alpha_type constexpr y1 = 3.978578321256066662987441200793155884L;
				alpha_type constexpr y2 = 0.859750287076040352418274981236440341L;
				alpha_type constexpr y3 = 0.022706582386768037080226242510354953L;
				auto const x = termial_f(w, x0, x1, x2, x3);
				auto const y = termial_f(w, y0, y1, y2, y3);
				return complexion_f(x, y);
			}
			XTAL_0IF (N_lim == 0x7) {// 8:7 D[...4]@0 && D[...5]@¼
				alpha_type constexpr x0 = one;
				alpha_type constexpr x1 = 4.606547402370024772724248519073301565L;
				alpha_type constexpr x2 = 2.500953736670776238483254500604284293L;
				alpha_type constexpr x3 = 0.300243870518009486091801923099845388L;
				alpha_type constexpr x4 = 0.004749448374844164315053279863807317L;
				alpha_type constexpr y0 = std::numbers::pi_v<alpha_type>;
				alpha_type constexpr y1 = 4.136469917598875065501883380414805209L;
				alpha_type constexpr y2 = 1.048974757862995076084471175276645114L;
				alpha_type constexpr y3 = 0.054095324243024101847331195900039153L;
				auto const x = termial_f(w, x0, x1, x2, x3, x4);
				auto const y = termial_f(w, y0, y1, y2, y3);
				return complexion_f(x, y);
			}
			XTAL_0IF_(else) {
				auto const [u, _u] = roots_f<2>(w);
				auto const x = cos(u*_fix::patio_1);
				auto const y = sin(u*_fix::patio_1)*_u;
				return complexion_f(x, y);
			}
		}

	};
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct impunity<M_ism,-2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(long,static)
		XTAL_LET function(simplex_field_q auto o)
		noexcept -> auto
		{
			using    _fix = bond::fixture<decltype(o)>;
			XTAL_LET _dn = one/_fix::patio_1;
			
			using U_aphex = typename _fix::aphex_type;
			using U_alpha = typename _fix::alpha_type;

			auto const w = objective_f(o);

			XTAL_IF0
			XTAL_0IF (0 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.18114206708921260e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 0.50387678776821820e-0;

				auto const x = termial_f(w, x0, x1);
				auto const y = termial_f(w, y0, y1);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (1 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.41502448325181734e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 0.74755993126133540e-0;
				XTAL_LET y2 =     (U_alpha) 0.05410519758331759e-0;
				
				auto const x = termial_f(w, x0, x1);
				auto const y = termial_f(w, y0, y1, y2);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (2 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.65051629987764580e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.03944474710871034e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 0.98379496428721650e-0;
				XTAL_LET y2 =     (U_alpha) 0.16793026979785060e-0;
				
				auto const x = termial_f(w, x0, x1, x2);
				auto const y = termial_f(w, y0, y1, y2);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (3 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.88375887183199590e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.13156850890207750e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.21708869440556770e-0;
				XTAL_LET y2 =     (U_alpha) 0.33731712522645585e-0;
				XTAL_LET y3 =     (U_alpha) 0.11588697106135937e-1;
				
				auto const x = termial_f(w, x0, x1, x2);
				auto const y = termial_f(w, y0, y1, y2, y3);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (4 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.10086520358652850e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.26813142431485410e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.07831573963866360e-1;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.43419816238477750e-0;
				XTAL_LET y2 =     (U_alpha) 0.54620442139995930e-0;
				XTAL_LET y3 =     (U_alpha) 0.45869073871868164e-1;
				
				auto const x = termial_f(w, x0, x1, x2, x3);
				auto const y = termial_f(w, y0, y1, y2, y3);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (5 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.33357813056109360e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.46539035057169720e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.35261607756788310e-1;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.66691144301400500e-0;
				XTAL_LET y2 =     (U_alpha) 0.82102801941969360e-0;
				XTAL_LET y3 =     (U_alpha) 1.18407591369199770e-1;
				XTAL_LET y4 =     (U_alpha) 0.23067742495692903e-2;
				
				auto const x = termial_f(w, x0, x1, x2, x3);
				auto const y = termial_f(w, y0, y1, y2, y3, y4);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (6 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.54923781439012060e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.69848490582322830e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.90506751518911740e-1;
				XTAL_LET x4 = _dn*(U_alpha) 0.15530589255192263e-2;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.88257114554573410e-0;
				XTAL_LET y2 =     (U_alpha) 1.12600868501140260e-0;
				XTAL_LET y3 =     (U_alpha) 2.32185184843683330e-1;
				XTAL_LET y4 =     (U_alpha) 1.15781734483417240e-2;

				auto const x = termial_f(w, x0, x1, x2, x3, x4);
				auto const y = termial_f(w, y0, y1, y2, y3, y4);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (7 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.78189081352097830e-0;
				XTAL_LET x2 = _dn*(U_alpha) 1.00031477243997680e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.19176306039358512e-0;
				XTAL_LET x4 = _dn*(U_alpha) 0.88074505355593190e-2;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 2.11522414674357060e-0;
				XTAL_LET y2 =     (U_alpha) 1.50538949216831600e-0;
				XTAL_LET y3 =     (U_alpha) 0.41337181267810030e-0;
				XTAL_LET y4 =     (U_alpha) 0.36584362720576720e-1;
				XTAL_LET y5 =     (U_alpha) 0.45821007587559800e-3;

				auto const x = termial_f(w, x0, x1, x2, x3, x4);
				auto const y = termial_f(w, y0, y1, y2, y3, y4, y5);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF_(else) {
				auto [d, q] = roots_f<2, 1>(o);
				(void) cut_t<XTAL_VAL_(-_fix::maxilon_1)>::edit(d);
				(void) cut_t<XTAL_VAL_(-_fix::maxilon_1)>::edit(q);
				XTAL_IF0
				XTAL_0IF (M_ism == -1) {return _dn*q*atan (d);}
				XTAL_0IF (M_ism == -2) {return _dn*q*atanh(d);}
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
