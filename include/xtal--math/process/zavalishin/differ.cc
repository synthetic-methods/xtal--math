#pragma once
#include "./any.cc"





#include "./differ.hh"
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _differ
{///////////////////////////////////////////////////////////////////////////////

struct reparameterized
{
	template <class S>
	class subtype : public S
	{
		using S_ = S;

	public:
		using S_::S_;

		using order_attribute = occur::inferred_t<union ORDER,
			bond::seek_s<1 + S_::data_type::size()>>;

	};
};


}///////////////////////////////////////////////////////////////////////////////

/**/
TAG_("differ")
{
	using U_fit = bond::fit<>;
	
	using          U_pole = typename U_fit::alpha_type;
	auto constexpr N_pole = fixed<U_pole[1]>::extent();

	//\
	using D1     = differ<_differ::reparameterized, U_pole[1]>;
	using D1     = differ<U_pole[1]>;
	using D2     = differ<U_pole[2]>;
	using D1_etc = occur::codex_t<D1>;
	using D2_etc = occur::codex_t<D2>;
	using D1_prx = process::confined_t<void
	,	typename D1_etc::template dispatch<>
	,	D1
	>;
	using D2_prx = process::confined_t<void
	,	typename D2_etc::template dispatch<>
	,	D2
	>;

	using U0_order = occur::inferred_t<union ORDER, bond::seek_s<1 + N_pole>>;
	using U1_order = typename D1_etc::order_attribute;
	using U2_order = typename D2_etc::order_attribute;
//	static_assert(same_q<U1_order, typename D1_prx::order_attribute>);
//	static_assert(same_q<U0_order, typename D1_prx::order_attribute>);

	TRY_("differentiation <N_sub=0>")
	{
		D1_prx d1{};
		D2_prx d2{};
		//\
		d1 <<= typename D1_prx::order_attribute{1};
		d1 <<= typename D1_etc::order_attribute{1};
		d2 <<= typename D2_etc::order_attribute{1};

		auto x = U_fit::haplo_f(7);
		auto y = U_fit::alpha_f(0);

		y = d1.template method<1, 0>(1*x); TRUE_(check_f<- 1>(x, y));
		y = d1.template method<1, 0>(2*x); TRUE_(check_f<- 1>(x, y));
		y = d1.template method<1, 0>(3*x); TRUE_(check_f<- 1>(x, y));
		y = d1.template method<1, 0>(4*x); TRUE_(check_f<- 1>(x, y));
		y = d1.template method<1, 0>(5*x); TRUE_(check_f<- 1>(x, y));
		y = d1.template method<1, 0>(6*x); TRUE_(check_f<- 1>(x, y));
		y = d1.template method<1, 0>(7*x); TRUE_(check_f<- 1>(x, y));
		y = d1.template method<1, 0>(8*x); TRUE_(check_f<- 1>(x, y));

		y = d2.template method<1, 0>(1*x); TRUE_(check_f<-44>(x, y));
		y = d2.template method<1, 0>(2*x); TRUE_(check_f<-43>(x, y));
		y = d2.template method<1, 0>(3*x); TRUE_(check_f<-42>(x, y));
		y = d2.template method<1, 0>(4*x); TRUE_(check_f<-40>(x, y));
		y = d2.template method<1, 0>(5*x); TRUE_(check_f<-39>(x, y));
		y = d2.template method<1, 0>(6*x); TRUE_(check_f<-38>(x, y));
		y = d2.template method<1, 0>(7*x); TRUE_(check_f<-37>(x, y));
		y = d2.template method<1, 0>(8*x); TRUE_(check_f<-35>(x, y));

	}
	TRY_("differentiation <N_sub=1>")
	{
		D1_prx d1{};
		D2_prx d2{};
		//\
		d1 <<= typename D1_prx::order_attribute{1};
		d1 <<= typename D1_etc::order_attribute{1};
		d2 <<= typename D2_etc::order_attribute{1};

		auto x = U_fit::haplo_f(7);
		auto y = U_fit::alpha_f(0);

		y = d1.template method<1, 1>(1*x); TRUE_(check_f<-44>(x, y));
		y = d1.template method<1, 1>(2*x); TRUE_(check_f<-43>(x, y));
		y = d1.template method<1, 1>(3*x); TRUE_(check_f<-42>(x, y));
		y = d1.template method<1, 1>(4*x); TRUE_(check_f<-40>(x, y));
		y = d1.template method<1, 1>(5*x); TRUE_(check_f<-39>(x, y));
		y = d1.template method<1, 1>(6*x); TRUE_(check_f<-38>(x, y));
		y = d1.template method<1, 1>(7*x); TRUE_(check_f<-37>(x, y));
		y = d1.template method<1, 1>(8*x); TRUE_(check_f<-35>(x, y));

		y = d2.template method<1, 1>(1*x); TRUE_(check_f<-45>(x, y));
		y = d2.template method<1, 1>(2*x); TRUE_(check_f<-44>(x, y));
		y = d2.template method<1, 1>(3*x); TRUE_(check_f<-43>(x, y));
		y = d2.template method<1, 1>(4*x); TRUE_(check_f<-40>(x, y));
		y = d2.template method<1, 1>(5*x); TRUE_(check_f<-41>(x, y));
		y = d2.template method<1, 1>(6*x); TRUE_(check_f<-41>(x, y));
		y = d2.template method<1, 1>(7*x); TRUE_(check_f<-40>(x, y));
		y = d2.template method<1, 1>(8*x); TRUE_(check_f<-39>(x, y));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
