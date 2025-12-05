#pragma once
#include "./any.cc"





#include "./arc.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("arc")
{
	using _fit = bond::fit<>;

	using U_sigma = typename _fit::sigma_type;
	using U_delta = typename _fit::delta_type;
	using U_alpha = typename _fit::alpha_type;
	using U_aphex = typename _fit::aphex_type;
	static constexpr U_alpha egg =  1.23456789;
	static constexpr U_alpha ten = 10;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	/**/
	TRY_("arc<-0> evaluation")
	{
		int  constexpr M_lim = -0;
		int  constexpr M_car =  0;
		auto constexpr fN = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f<~0   >);
		auto constexpr f7 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 2>);
		auto constexpr f6 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 1>);
		auto constexpr f5 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 2>);
		auto constexpr f4 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 1>);
		auto constexpr f3 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 2>);
		auto constexpr f2 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 1>);
		auto constexpr f1 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 2>);
		auto constexpr f0 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 1>);

		U_alpha t{0.123};
		U_alpha u{0.456};
		U_alpha v{0.789};
		U_aphex w{u, v};

		TRUE_(check_f<-26>(fN(u), f3(u)));
		TRUE_(check_f<-31>(fN(u), f2(u)));
		TRUE_(check_f<-36>(fN(u), f1(u)));
		TRUE_(check_f<-42>(fN(u), f0(u)));

		TRUE_(check_f<-32>(fN(v), f3(v)));
		TRUE_(check_f<-36>(fN(v), f2(v)));
		TRUE_(check_f<-40>(fN(v), f1(v)));
		TRUE_(check_f<-44>(fN(v), f0(v)));

		TRUE_(check_f<-41>(fN(w), f3(w)));
		TRUE_(check_f<-43>(fN(w), f2(w)));
		TRUE_(check_f<-45>(fN(w), f1(w)));
		TRUE_(check_f<-47>(fN(w), f0(w)));

		TRUE_(check_f<-24>(fN(t/w), f3(t, w)));
		TRUE_(check_f<-29>(fN(t/w), f2(t, w)));
		TRUE_(check_f<-33>(fN(t/w), f1(t, w)));
		TRUE_(check_f<-38>(fN(t/w), f0(t, w)));

	}
	/***/
	/**/
	TRY_("arc<-0> evaluation w/ argument restriction")
	{
		int  constexpr M_lim = -0;
		int  constexpr M_car =  1;
		auto constexpr fN = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f<~0   >);
		auto constexpr f7 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 2>);
		auto constexpr f6 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 1>);
		auto constexpr f5 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 2>);
		auto constexpr f4 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 1>);
		auto constexpr f3 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 2>);
		auto constexpr f2 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 1>);
		auto constexpr f1 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 2>);
		auto constexpr f0 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 1>);

		U_alpha t{2.123};
		U_alpha u{5.456};
		U_alpha v{8.789};
		U_aphex w{u, v};

		TRUE_(check_f<- 1>(fN(u), f7(u)));
		TRUE_(check_f<- 6>(fN(u), f6(u)));
		TRUE_(check_f<-12>(fN(u), f5(u)));
		TRUE_(check_f<-17>(fN(u), f4(u)));
		TRUE_(check_f<-21>(fN(u), f3(u)));
		TRUE_(check_f<-28>(fN(u), f2(u)));
		TRUE_(check_f<-34>(fN(u), f1(u)));
		TRUE_(check_f<-39>(fN(u), f0(u)));
	
		TRUE_(check_f<- 1>(fN(v), f7(v)));
		TRUE_(check_f<- 1>(fN(v), f6(v)));
		TRUE_(check_f<- 9>(fN(v), f5(v)));
		TRUE_(check_f<-16>(fN(v), f4(v)));
		TRUE_(check_f<-22>(fN(v), f3(v)));
		TRUE_(check_f<-27>(fN(v), f2(v)));
		TRUE_(check_f<-32>(fN(v), f1(v)));
		TRUE_(check_f<-37>(fN(v), f0(v)));
		
	}
	/***/

	/**/
	TRY_("arc<-1> evaluation")
	{
		int  constexpr M_lim = -1;
		int  constexpr M_car =  0;
		auto constexpr fN = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f<~0   >);
		auto constexpr f7 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 2>);
		auto constexpr f6 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 1>);
		auto constexpr f5 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 2>);
		auto constexpr f4 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 1>);
		auto constexpr f3 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 2>);
		auto constexpr f2 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 1>);
		auto constexpr f1 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 2>);
		auto constexpr f0 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 1>);

		U_alpha t{0.123};
		U_alpha u{0.456};
		U_alpha v{0.789};
		U_aphex w{u, v};

		TRUE_(check_f<-19>(fN(u), f7(u)));
		TRUE_(check_f<-21>(fN(u), f6(u)));
		TRUE_(check_f<-24>(fN(u), f5(u)));
		TRUE_(check_f<-28>(fN(u), f4(u)));
		TRUE_(check_f<-31>(fN(u), f3(u)));
		TRUE_(check_f<-36>(fN(u), f2(u)));
		TRUE_(check_f<-19>(fN(u), f1(u)));
		TRUE_(check_f<-43>(fN(u), f0(u)));

		TRUE_(check_f<- 8>(fN(v), f7(v)));
		TRUE_(check_f<- 9>(fN(v), f6(v)));
		TRUE_(check_f<-14>(fN(v), f5(v)));
		TRUE_(check_f<-19>(fN(v), f4(v)));
		TRUE_(check_f<-24>(fN(v), f3(v)));
		TRUE_(check_f<-27>(fN(v), f2(v)));
		TRUE_(check_f<- 8>(fN(v), f1(v)));
		TRUE_(check_f<-38>(fN(v), f0(v)));

	}
	/**/
	TRY_("arc<-1> evaluation w/ argument restriction")
	{
		int  constexpr M_lim = -1;
		int  constexpr M_car =  1;
		auto constexpr fN = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f<~0   >);
		auto constexpr f7 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 2>);
		auto constexpr f6 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 3, 1>);
		auto constexpr f5 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 2>);
		auto constexpr f4 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 2, 1>);
		auto constexpr f3 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 2>);
		auto constexpr f2 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 1, 1>);
		auto constexpr f1 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 2>);
		auto constexpr f0 = [] XTAL_1FN_(call) (arc_t<M_lim, M_car>::template method_f< 0, 1>);

		U_alpha t{0.123};
		U_alpha u{0.456};
		U_alpha v{0.789};
		U_aphex w{u, v};

		TRUE_(check_f<-23>(fN(u), f3(u)));
		TRUE_(check_f<-28>(fN(u), f2(u)));
		TRUE_(check_f<- 8>(fN(u), f1(u)));
		TRUE_(check_f<-37>(fN(u), f0(u)));

		TRUE_(check_f<-20>(fN(v), f3(v)));
		TRUE_(check_f<-25>(fN(v), f2(v)));
		TRUE_(check_f<- 3>(fN(v), f1(v)));
		TRUE_(check_f<-33>(fN(v), f0(v)));

	}
	/***/
};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
