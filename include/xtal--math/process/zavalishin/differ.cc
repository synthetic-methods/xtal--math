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
	public:
		using S::S;

		using typename S::pole_size;
		using order_type = occur::inferred_t<union ORDER, bond::seek_s<1 + pole_size{}()>>;

	};
};


}///////////////////////////////////////////////////////////////////////////////

/**/
TAG_("differ")
{
	using          U_pole = typename bond::fit<>::alpha_type;
	auto constexpr N_pole = fixed_n<U_pole[1]>;

	//\
	using D_def = differ<_differ::reparameterized, U_pole[1]>;
	using D_def = differ<U_pole[1]>;
	using D_etc = process::traits_t<D_def>;
	using D_prx = process::confined_t<void
	,	typename D_etc::template dispatch<>
	,	D_def
	>;

	using V_order = occur::inferred_t<union ORDER, bond::seek_s<1 + N_pole>>;
	using U_order = typename D_etc::order_type;
//	static_assert(same_q<U_order, typename D_prx::order_type>);
//	static_assert(same_q<V_order, typename D_prx::order_type>);

	TRY_("instantiation")
	{
		D_prx diff{};
		//\
		diff <<= typename D_prx::order_type{1};
		diff <<= typename D_etc::order_type{1};

		TRUE_(check_f<19>(diff(1), 1));
		TRUE_(check_f<19>(diff(2), 1));
		TRUE_(check_f<19>(diff(4), 2));
		TRUE_(check_f<19>(diff(8), 4));

	}
//	TRY_("scaling")
//	{
//		confined_t<dilate<[] XTAL_1FN_(to) (bond::fit<>::patio_2)>, D_prx> diff{};
//		diff <<= typename D_prx::order_type{1};
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
