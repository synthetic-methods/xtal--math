#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_rot=0, int M_car=0>
struct  arc;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Rational approximation for the arctangent function `ArcTan[#]`.
*/
template <>
struct arc<-0, 0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim, N_par>(XTAL_REF_(o), unstruct_t<XTAL_ALL_(o)>{one});
		}
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&y, auto &&x)
		noexcept -> decltype(auto)
		{
			using          U_fit = bond::fit<decltype(y), decltype(x)>;
			auto constexpr  _1pi = one/U_fit::patio_1;
			auto constexpr  _2pi = one/U_fit::patio_2;
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)   {return mythos_f<N_lim, N_par>(XTAL_REF_(y), XTAL_REF_(x));}
			XTAL_0IF_(consteval)    {return mythos_f<   ~0, N_par>(XTAL_REF_(y), XTAL_REF_(x));}
	#if   XTAL_SYS_(builtin)
			XTAL_0IF (real_variable_q<decltype(y)> and real_variable_q<decltype(x)>)
			                        {return  __builtin_atan(XTAL_REF_(y)/XTAL_REF_(x))*_1pi;}
	#endif
			XTAL_0IF_(else)         {return            atan(XTAL_REF_(y)/XTAL_REF_(x))*_1pi;}
		}

	protected:
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		mythos_f(auto &&y, auto &&x)
		noexcept -> decltype(auto)
		{
			auto constexpr I_par = 1 + (~N_par&1);
			auto constexpr I_lim = below_v<4, (unsigned) N_lim>;
			auto     const yy =   square_f(y);
			auto     const xx =   square_f(x);
			using           X =  XTAL_ALL_(x);
			using           Y =  XTAL_ALL_(y);
			using           K = unstruct_t<Y>;
			K constexpr   _pi = 0.318309886183790671537767526745029L;// 1/Pi
			auto const   y_pi = XTAL_REF_(y)*_pi;
			XTAL_IF0
		//	N_par == 1
			XTAL_0IF (I_lim == 0 and I_par == 1) {
				K constexpr X0 = 0.273239544735162686151070106980115L;// 4/Pi - 1
				auto  const dn = term_f(xx, yy, X0);
				return           root_f<-1>(dn) * (XTAL_REF_(x)*XTAL_MOV_(y_pi));
			}
			XTAL_0IF (I_lim == 1 and I_par == 1) {
				K constexpr X0 = 0.092632089554137291566353742342042L, Y0 = 0.455708497276173034068352626103896L;
				K constexpr X1 = 0.696330944385271600741835258964320L;
				auto  const dn = term_f(xx, yy, X0)*
				                 term_f(xx, yy, X1);
				auto  const up = term_f(xx, yy, Y0);
				return      up * root_f<-1>(dn) * (XTAL_REF_(x)*XTAL_MOV_(y_pi));
			}
			XTAL_0IF (I_lim == 2 and I_par == 1) {
				K constexpr X0 = 0.045336879390751627376445716562634L, Y0 = 0.724293136180053939407773098279180L;
				K constexpr X1 = 0.382148801370432439320822739969860L, Y1 = 0.212898303369872490915300981673779L;
				K constexpr X2 = 0.843039031089519011457888726457800L;
				auto  const dn = term_f(xx, yy, X0)*
				                 term_f(xx, yy, X1)*
				                 term_f(xx, yy, X2);
				auto  const up = term_f(xx, yy, Y0)*
				                 term_f(xx, yy, Y1);
				return      up * root_f<-1>(dn) * (XTAL_REF_(x)*XTAL_MOV_(y_pi));
			}
			XTAL_0IF (I_lim == 3 and I_par == 1) {
				K constexpr X0 = 0.026732055391081833191956738145699L, Y0 = 0.839759054106618847678222417078600L;
				K constexpr X1 = 0.232643129601723913230550038619120L, Y1 = 0.452699143092020462022097547694890L;
				K constexpr X2 = 0.582303363755852685124094797824660L, Y2 = 0.121330001415257447615383301822009L;
				K constexpr X3 = 0.905442983158669863808040841639100L;
				auto  const dn = term_f(xx, yy, X0)*
				                 term_f(xx, yy, X1)*
				                 term_f(xx, yy, X2)*
				                 term_f(xx, yy, X3);
				auto  const up = term_f(xx, yy, Y0)*
				                 term_f(xx, yy, Y1)*
				                 term_f(xx, yy, Y2);
				return      up * root_f<-1>(dn) * (XTAL_REF_(x)*XTAL_MOV_(y_pi));
			}
		//	N_par == 2
			XTAL_0IF (I_lim == 0 and I_par == 2) {
				K constexpr X0 = 0.541799971604698858679141097772248L, Y0 = 0.210926866024568438254555477047590L;
				auto  const dn = term_f(xx, yy, X0);
				auto  const up = term_f(xx, yy, Y0);
				return      up * root_f<-1>(dn*XTAL_REF_(x)) * (XTAL_MOV_(y_pi));
			}
			XTAL_0IF (I_lim == 1 and I_par == 2) {
				K constexpr X0 = 0.786808235797912385470876726654560L, Y0 = 0.619519329848094497184408154973100L;
				K constexpr X1 = 0.244005542165143401309672324805773L, Y1 = 0.077963376810944219129679487901419L;
				auto  const dn = term_f(xx, yy, X0)*
				                 term_f(xx, yy, X1);
				auto  const up = term_f(xx, yy, Y0)*
				                 term_f(xx, yy, Y1);
				return      up * root_f<-1>(dn*XTAL_REF_(x)) * (XTAL_MOV_(y_pi));
			}
			XTAL_0IF (I_lim == 2 and I_par == 2) {
				K constexpr X0 = 0.134918588688084249936941384877666L, Y0 = 0.793008417546886251007533332037000L;
				K constexpr X1 = 0.494230041597349895472205531005920L, Y1 = 0.342786492530276980497987158197830L;
				K constexpr X2 = 0.879993343477236904116537635581440L, Y2 = 0.040013731934545967039340988888031L;
				auto  const dn = term_f(xx, yy, X0)*
				                 term_f(xx, yy, X1)*
				                 term_f(xx, yy, X2);
				auto  const up = term_f(xx, yy, Y0)*
				                 term_f(xx, yy, Y1)*
				                 term_f(xx, yy, Y2);
				return      up * root_f<-1>(dn*XTAL_REF_(x)) * (XTAL_MOV_(y_pi));
			}
			//\
			XTAL_0IF (I_lim == 3 and I_par == 2) {
			XTAL_0IF_(else) {
				K constexpr X0 = 0.085067308035404941721715359954706L, Y0 = 0.024250681053273100097897318849783L;
				K constexpr X1 = 0.324760285686760780911158131455339L, Y1 = 0.212487646203152306639435622436891L;
				K constexpr X2 = 0.651209859798724551097201190480120L, Y2 = 0.541909515304353095181096867463551L;
				K constexpr X3 = 0.923658437079931678056766743691228L, Y3 = 0.872714714707695932850422917475650L;
				auto  const dn = term_f(xx, yy, X0)*
				                 term_f(xx, yy, X1)*
				                 term_f(xx, yy, X2)*
				                 term_f(xx, yy, X3);
				auto  const up = term_f(xx, yy, Y0)*
				                 term_f(xx, yy, Y1)*
				                 term_f(xx, yy, Y2)*
				                 term_f(xx, yy, Y3);
				return      up * root_f<-1>(dn*XTAL_REF_(x)) * (XTAL_MOV_(y_pi));
			}
		}

	};
};
template <>
struct arc<-0, 1>
{
	using superkind = arc<-0, 0>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim, N_par>(XTAL_REF_(o), unstruct_t<XTAL_ALL_(o)>{one});
		}
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&y, auto &&x)
		noexcept -> decltype(auto)
		{
			using         U_fit = bond::fit<decltype(y), decltype(x)>;
			auto constexpr _1pi = one/U_fit::patio_1;
			auto constexpr _2pi = one/U_fit::patio_2;
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)   {return mythos_f<N_lim, N_par>(XTAL_REF_(y), XTAL_REF_(x));}
			XTAL_0IF_(consteval)    {return mythos_f<   ~0, N_par>(XTAL_REF_(y), XTAL_REF_(x));}
	#if   XTAL_SYS_(builtin)
			XTAL_0IF (real_variable_q<decltype(y)> and real_variable_q<decltype(x)>)
			                        {return  __builtin_atan(XTAL_REF_(y)/XTAL_REF_(x))*_1pi;}
	#endif
			XTAL_0IF_(else)         {return            atan(XTAL_REF_(y)/XTAL_REF_(x))*_1pi;}
		}

	protected:
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,set)
		mythos_f(auto &&v, auto &&u)
		noexcept -> auto
		{
			static_assert(real_variable_q<decltype(v)>);
			static_assert(real_variable_q<decltype(u)>);
			using L = bond::fit<decltype(v), decltype(u)>;
			using U_aphex = typename L::aphex_type;
			using U_alpha = typename L::alpha_type;
			using W_alpha = atom::couple_t<U_alpha[2]>;

			auto u_abs = u, u_sgn = decompose_t<signed>::method_e(u_abs);
			auto v_abs = v, v_sgn = decompose_t<signed>::method_e(v_abs);// v_sgn *= *L::haplo_1;

			W_alpha co{v_abs < u_abs, _std::in_place};
			W_alpha up{v, u_abs}; up *= co;
			W_alpha dn{u_abs,-v}; dn *= co;

			auto const &[co_0, co_1]  = co;
			auto const u_flp = L::haplo_1 - L::haplo_1*co_0*u_sgn;

			return term_f(u_flp*v_sgn, u_sgn, S_::template method_f<N_lim, N_par>(up.sum(), dn.sum()));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Rational approximation for the logarithm function.
*/
template <>
struct arc<-1, 0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			using L = bond::fit<decltype(o)>;
			auto constexpr  _1pi = one/L::patio_f(-1);
			auto constexpr  _2pi = one/L::patio_f(-2);
			auto constexpr  _4pi = one/L::patio_f(-4);
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)                   {return mythos_f<N_lim, N_par>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)                    {return mythos_f<   ~0, N_par>(XTAL_REF_(o));}
	#if   XTAL_SYS_(builtin)
			XTAL_0IF (real_variable_q<decltype(o)>) {return     _2pi*__builtin_log(XTAL_REF_(o));}
	#endif
			XTAL_0IF_(else)                         {return     _2pi*          log(XTAL_REF_(o));}
		}
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&y, auto &&x)
		noexcept -> decltype(auto)
		{
			return half*method_f<N_lim, N_par>(square_f(XTAL_REF_(y), XTAL_REF_(x)));
		}

	protected:
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		mythos_f(auto &&u)
		noexcept -> decltype(auto)
		{
			auto constexpr I_par = 1 + (~N_par&1);
			auto constexpr I_lim = below_v<4, (unsigned) N_lim>;
			using K = unstruct_t<XTAL_ALL_(u)>;
			using L = bond::fit<K>;
			auto constexpr  _1pi = one/L::patio_f(-1);
			auto constexpr  _2pi = one/L::patio_f(-2);
			auto constexpr  _4pi = one/L::patio_f(-4);
			auto  const w    = square_f(u);
			auto  const w_up = w + one, w_dn = w - one;
			auto  const u_up = u + one, u_dn = u - one;
			XTAL_IF0
		//	N_par == 1
			XTAL_0IF (I_lim == 0 and I_par == 1) {
				K constexpr    X0 = 0.2732395447351626861510701069801150L;// 4/Pi - 1
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				auto const  dn    = term_f(u, w_up, dn_X0/up_X0/2);
				auto const  up    =                  _2pi/up_X0;
				return    w_dn*up * root_f<-1>(dn);
			}
			XTAL_0IF (I_lim == 1 and I_par == 1) {
				K constexpr    X0 = 0.092632089554137291566353742342042L, Y0 = 0.455708497276173034068352626103896L;
				K constexpr    X1 = 0.696330944385271600741835258964320L;
				K constexpr dn_Y0 = 1 - Y0, up_Y0 = 1 + Y0;
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				K constexpr dn_X1 = 1 - X1, up_X1 = 1 + X1;
				K constexpr dn_co =         up_X0*
				                            up_X1;
				K constexpr up_co =         up_Y0*_2pi;
				auto  const    dn = term_f(u, w_up, dn_X0/up_X0/2)*
				                    term_f(u, w_up, dn_X1/up_X1/2);
				auto  const    up = term_f(u, w_up, dn_Y0/up_Y0/2);
				return    w_dn*up * root_f<-1>(dn) * (up_co/dn_co);
			}
			XTAL_0IF (I_lim == 2 and I_par == 1) {
				K constexpr    X0 = 0.045336879390751627376445716562634L, Y0 = 0.724293136180053939407773098279180L;
				K constexpr    X1 = 0.382148801370432439320822739969860L, Y1 = 0.212898303369872490915300981673779L;
				K constexpr    X2 = 0.843039031089519011457888726457800L;
				K constexpr dn_Y0 = 1 - Y0, up_Y0 = 1 + Y0;
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				K constexpr dn_Y1 = 1 - Y1, up_Y1 = 1 + Y1;
				K constexpr dn_X1 = 1 - X1, up_X1 = 1 + X1;
				K constexpr dn_X2 = 1 - X2, up_X2 = 1 + X2;
				K constexpr dn_co =         up_X0*
				                            up_X1*
				                            up_X2;
				K constexpr up_co =         up_Y0*
				                            up_Y1*_2pi;
				auto  const    dn = term_f(u, w_up, dn_X0/up_X0/2)*
				                    term_f(u, w_up, dn_X1/up_X1/2)*
				                    term_f(u, w_up, dn_X2/up_X2/2);
				auto  const    up = term_f(u, w_up, dn_Y0/up_Y0/2)*
				                    term_f(u, w_up, dn_Y1/up_Y1/2);
				return    w_dn*up * root_f<-1>(dn) * (up_co/dn_co);
			}
			XTAL_0IF (I_lim == 3 and I_par == 1) {
				K constexpr    X0 = 0.026732055391081833191956738145699L, Y0 = 0.839759054106618847678222417078600L;
				K constexpr    X1 = 0.232643129601723913230550038619120L, Y1 = 0.452699143092020462022097547694890L;
				K constexpr    X2 = 0.582303363755852685124094797824660L, Y2 = 0.121330001415257447615383301822009L;
				K constexpr    X3 = 0.905442983158669863808040841639100L;
				K constexpr dn_Y0 = 1 - Y0, up_Y0 = 1 + Y0;
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				K constexpr dn_Y1 = 1 - Y1, up_Y1 = 1 + Y1;
				K constexpr dn_X1 = 1 - X1, up_X1 = 1 + X1;
				K constexpr dn_X2 = 1 - X2, up_X2 = 1 + X2;
				K constexpr dn_Y2 = 1 - Y2, up_Y2 = 1 + Y2;
				K constexpr dn_X3 = 1 - X3, up_X3 = 1 + X3;
				K constexpr dn_co =         up_X0*
				                            up_X1*
				                            up_X2*
				                            up_X3;
				K constexpr up_co =         up_Y0*
				                            up_Y1*
				                            up_Y2*_2pi;
				auto  const    dn = term_f(u, w_up, dn_X0/up_X0/2)*
				                    term_f(u, w_up, dn_X1/up_X1/2)*
				                    term_f(u, w_up, dn_X2/up_X2/2)*
				                    term_f(u, w_up, dn_X3/up_X3/2);
				auto  const    up = term_f(u, w_up, dn_Y0/up_Y0/2)*
				                    term_f(u, w_up, dn_Y1/up_Y1/2)*
				                    term_f(u, w_up, dn_Y2/up_Y2/2);
				return    w_dn*up * root_f<-1>(dn) * (up_co/dn_co);
			}
		//	N_par == 2
			XTAL_0IF (I_lim == 0 and I_par == 1) {
				K constexpr    X0 = 0.541799971604698858679141097772248L, Y0 = 0.210926866024568438254555477047590L;
				K constexpr dn_Y0 = 1 - Y0, up_Y0 = 1 + Y0;
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				K constexpr dn_co =         up_X0;
				K constexpr up_co =         up_Y0*_1pi;
				auto const     up = term_f(u, w_up, dn_Y0/up_Y0/2);
				auto const     dn = term_f(u, w_up, dn_X0/up_X0/2);
				return    u_dn*up * root_f<-1>(u_up*dn) * (up_co/dn_co);
			}
			XTAL_0IF (I_lim == 1 and I_par == 2) {
				K constexpr    X0 = 0.786808235797912385470876726654560L, Y0 = 0.619519329848094497184408154973100L;
				K constexpr    X1 = 0.244005542165143401309672324805773L, Y1 = 0.077963376810944219129679487901419L;
				K constexpr dn_Y0 = 1 - Y0, up_Y0 = 1 + Y0;
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				K constexpr dn_Y1 = 1 - Y1, up_Y1 = 1 + Y1;
				K constexpr dn_X1 = 1 - X1, up_X1 = 1 + X1;
				K constexpr dn_co =         up_X0*
				                            up_X1;
				K constexpr up_co =         up_Y0*
				                            up_Y1*_1pi;
				auto const     dn = term_f(u, w_up, dn_X0/up_X0/2)*
				                    term_f(u, w_up, dn_X1/up_X1/2);
				auto const     up = term_f(u, w_up, dn_Y0/up_Y0/2)*
				                    term_f(u, w_up, dn_Y1/up_Y1/2);
				return    u_dn*up * root_f<-1>(u_up*dn) * (up_co/dn_co);
			}
			XTAL_0IF (I_lim == 2 and I_par == 2) {
				K constexpr    X0 = 0.134918588688084249936941384877666L, Y0 = 0.793008417546886251007533332037000L;
				K constexpr    X1 = 0.494230041597349895472205531005920L, Y1 = 0.342786492530276980497987158197830L;
				K constexpr    X2 = 0.879993343477236904116537635581440L, Y2 = 0.040013731934545967039340988888031L;
				K constexpr dn_Y0 = 1 - Y0, up_Y0 = 1 + Y0;
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				K constexpr dn_Y1 = 1 - Y1, up_Y1 = 1 + Y1;
				K constexpr dn_X1 = 1 - X1, up_X1 = 1 + X1;
				K constexpr dn_Y2 = 1 - Y2, up_Y2 = 1 + Y2;
				K constexpr dn_X2 = 1 - X2, up_X2 = 1 + X2;
				K constexpr dn_co =         up_X0*
				                            up_X1*
				                            up_X2;
				K constexpr up_co =         up_Y0*
				                            up_Y1*
				                            up_Y2*_1pi;
				auto const     dn = term_f(u, w_up, dn_X0/up_X0/2)*
				                    term_f(u, w_up, dn_X1/up_X1/2)*
				                    term_f(u, w_up, dn_X2/up_X2/2);
				auto const     up = term_f(u, w_up, dn_Y0/up_Y0/2)*
				                    term_f(u, w_up, dn_Y1/up_Y1/2)*
				                    term_f(u, w_up, dn_Y2/up_Y2/2);
				return    u_dn*up * root_f<-1>(u_up*dn) * (up_co/dn_co);
			}
			//\
			XTAL_0IF (I_lim == 3 and I_par == 2) {
			XTAL_0IF_(else) {
				K constexpr    X0 = 0.081864856415554748051056555790168L, Y0 = 0.023221333952416495390547391924276L;
				K constexpr    X1 = 0.315616185530778780488046972975596L, Y1 = 0.205078096516516402405519469750355L;
				K constexpr    X2 = 0.641627423211762055523958248087220L, Y2 = 0.530634436337003641943997201262050L;
				K constexpr    X3 = 0.920664457245484323471703998224308L, Y3 = 0.867505722294008916548053321941278L;
				K constexpr dn_Y0 = 1 - Y0, up_Y0 = 1 + Y0;
				K constexpr dn_X0 = 1 - X0, up_X0 = 1 + X0;
				K constexpr dn_Y1 = 1 - Y1, up_Y1 = 1 + Y1;
				K constexpr dn_X1 = 1 - X1, up_X1 = 1 + X1;
				K constexpr dn_Y2 = 1 - Y2, up_Y2 = 1 + Y2;
				K constexpr dn_X2 = 1 - X2, up_X2 = 1 + X2;
				K constexpr dn_Y3 = 1 - Y3, up_Y3 = 1 + Y3;
				K constexpr dn_X3 = 1 - X3, up_X3 = 1 + X3;
				K constexpr dn_co =         up_X0*
				                            up_X1*
				                            up_X2*
				                            up_X3;
				K constexpr up_co =         up_Y0*
				                            up_Y1*
				                            up_Y2*
				                            up_Y3*_1pi;
				auto const     dn = term_f(u, w_up, dn_X0/up_X0/2)*
				                    term_f(u, w_up, dn_X1/up_X1/2)*
				                    term_f(u, w_up, dn_X2/up_X2/2)*
				                    term_f(u, w_up, dn_X3/up_X3/2);
				auto const     up = term_f(u, w_up, dn_Y0/up_Y0/2)*
				                    term_f(u, w_up, dn_Y1/up_Y1/2)*
				                    term_f(u, w_up, dn_Y2/up_Y2/2)*
				                    term_f(u, w_up, dn_Y3/up_Y3/2);
				return    u_dn*up * root_f<-1>(u_up*dn) * (up_co/dn_co);
			}
		}

	};
};
template <>
struct arc<-1, 1>
{
	using superkind = arc<-1, 0>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> XTAL_ALL_(o)
		{
			using L = bond::fit<decltype(o)>;
			auto constexpr  _1pi = one/L::patio_f(-1);
			auto constexpr  _2pi = one/L::patio_f(-2);
			auto constexpr  _4pi = one/L::patio_f(-4);
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return methodology_f<N_lim, N_par>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)  {return methodology_f<   ~0, N_par>(XTAL_REF_(o));}
			XTAL_0IF_(else)       {return                    _2pi*log(XTAL_REF_(o));}
		}
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&y, auto &&x)
		noexcept -> decltype(auto)
		{
			return half*method_f<N_lim, N_par>(square_f(XTAL_REF_(y), XTAL_REF_(x)));
		}

	protected:
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,set)
		methodology_f(real_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using L = bond::fit<decltype(o)>;
			auto constexpr  _1pi = one/L::patio_f(-1);
			auto constexpr  _2pi = one/L::patio_f(-2);
			auto constexpr  _4pi = one/L::patio_f(-4);
			using U_alpha = typename L::alpha_type;
			using U_sigma = typename L::sigma_type;
			using U_delta = typename L::delta_type;

			U_alpha constexpr K_inv =      0.7071067811865475244008443621048490e0;//   1/Sqrt@2
			U_alpha constexpr K_log = _2pi*0.3465735902799726547086160607290883e0;// Log@Sqrt@2 // 2 Pi #&

		//	Log[m * 2^n] == Log[m] + Log[Sqrt@2]*n
		//	Log[m/Sqrt@2 * Sqrt@2^(1 + 2x)] == Log[m/Sqrt@2] + Log[Sqrt@2]*(1 + 2*n)*
			U_sigma m = _xtd::bit_cast<U_sigma>(o);
			U_delta n = m - L::unit.mask;
			m  &= L::fraction.mask;
			m  |= L::unit.mask;
			n >>= L::unit.shift - one;
			n  |=                 one;
			return term_f(S_::template method_f<N_lim, N_par>(K_inv*_xtd::bit_cast<U_alpha>(m)),
				K_log, static_cast<U_alpha>(n));
		}
		template <int N_lim=-1, int N_par=1>
		XTAL_DEF_(return,set)
		methodology_f(complex_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			auto const o_re = o.real();
			auto const o_im = o.imag();
			return complexion_f(method_f<N_lim, N_par>(o_im, o_re),
				confined_t<arc<0, 0>>::template method_f<N_lim, N_par>(o_im, o_re));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_rot=0, int M_car=0>
XTAL_TYP_(let) arc_t = process::confined_t<arc<M_rot, M_car>>;

template <int M_rot=0, int M_car=0, int ...Ns>
XTAL_DEF_(let) arc_f = [] XTAL_1FN_(call) (arc_t<M_rot, M_car>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
