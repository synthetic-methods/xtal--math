#pragma once
#include "./any.cc"





#include "./differ.hh"// testing...
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
			bond::seek_s<1 + S_::state_type::size()>>;

	};
};


}///////////////////////////////////////////////////////////////////////////////

/**/
TAG_("differ")
{
	using          U_pole = typename bond::fit<>::alpha_type;
	auto constexpr N_pole = fixed<U_pole[1]>::extent();

	//\
	using D_def = differ<_differ::reparameterized, U_pole[1]>;
	using D_def = differ<U_pole[1]>;
	using D_etc = occur::context_t<D_def>;
	using D_prx = process::confined_t<void
	,	typename D_etc::template dispatch<>
	,	D_def
	>;

	using V_order = occur::inferred_t<union ORDER, bond::seek_s<1 + N_pole>>;
	using U_order = typename D_etc::order_attribute;
//	static_assert(same_q<U_order, typename D_prx::order_attribute>);
//	static_assert(same_q<V_order, typename D_prx::order_attribute>);

	TRY_("instantiation")
	{
		D_prx diff{};
		//\
		diff <<= typename D_prx::order_attribute{1};
		diff <<= typename D_etc::order_attribute{1};

		TRUE_(check_f<19>(diff(1), 1));
		TRUE_(check_f<19>(diff(2), 1));
		TRUE_(check_f<19>(diff(4), 2));
		TRUE_(check_f<19>(diff(8), 4));

	}
//	TRY_("scaling")
//	{
//		confined_t<dilate<[] XTAL_1FN_(to) (bond::fit<>::patio_2)>, D_prx> diff{};
//		diff <<= typename D_prx::order_attribute{1};
//
//		TRUE_(check_f<19>(0.15915493667125702, diff(1)));
//		TRUE_(check_f<19>(0.15915493667125702, diff(2)));
//
//	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
