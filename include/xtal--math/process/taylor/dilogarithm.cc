#pragma once
#include "./any.cc"

#include "./tangent.hh"



#include "./dilogarithm.hh"
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("dilogarithm")
{
	using U_fit = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;
	static constexpr U_alpha egg = 0.123456789;
	static constexpr U_alpha ten = 10;

	using U_phi = atom::math::phason_t<U_alpha[2]>;

	/**/
	TRY_("evaluation")
	{
		TRUE_(check_f<-33>(dilogarithm_f<1, 0, 3>(U_aphex{0.77, 0.77}), U_aphex{ 0.5862403562859301, 1.0866295503986803}));
		TRUE_(check_f<-36>(dilogarithm_f<1, 0, 2>(U_aphex{0.77, 0.77}), U_aphex{ 0.5862403562859301, 1.0866295503986803}));
		TRUE_(check_f<-39>(dilogarithm_f<1, 0, 1>(U_aphex{0.77, 0.77}), U_aphex{ 0.5862403562859301, 1.0866295503986803}));
		TRUE_(check_f<-41>(dilogarithm_f<1, 0, 0>(U_aphex{0.77, 0.77}), U_aphex{ 0.5862403562859301, 1.0866295503986803}));
		
		TRUE_(check_f<-21>(dilogarithm_f<1, 0, 3>(U_aphex{0.22, 0.22}), U_aphex{ 0.21697134033239715, 0.2464548335260734}));
		TRUE_(check_f<-23>(dilogarithm_f<1, 0, 2>(U_aphex{0.22, 0.22}), U_aphex{ 0.21697134033239715, 0.2464548335260734}));
		TRUE_(check_f<-29>(dilogarithm_f<1, 0, 1>(U_aphex{0.22, 0.22}), U_aphex{ 0.21697134033239715, 0.2464548335260734}));
		TRUE_(check_f<-33>(dilogarithm_f<1, 0, 0>(U_aphex{0.22, 0.22}), U_aphex{ 0.21697134033239715, 0.2464548335260734}));
		
		TRUE_(check_f<-23>(dilogarithm_f<1, 0, 3>(U_aphex{0.77, 0.77}, U_aphex{0.11, 0.11}), U_aphex{-0.169842835980408100,-0.16886537033895920}));
		TRUE_(check_f<-23>(dilogarithm_f<1, 0, 2>(U_aphex{0.77, 0.77}, U_aphex{0.11, 0.11}), U_aphex{-0.169842835980408100,-0.16886537033895920}));
		TRUE_(check_f<-26>(dilogarithm_f<1, 0, 1>(U_aphex{0.77, 0.77}, U_aphex{0.11, 0.11}), U_aphex{-0.169842835980408100,-0.16886537033895920}));
		TRUE_(check_f<-42>(dilogarithm_f<1, 0, 0>(U_aphex{0.77, 0.77}, U_aphex{0.11, 0.11}), U_aphex{-0.169842835980408100,-0.16886537033895920}));

		TRUE_(check_f<-24>(dilogarithm_f<1, 0, 3>(U_aphex{0.22, 0.22}, U_aphex{0.11, 0.11}), U_aphex{-0.048411922454435596,-0.04838741284417837}));
		TRUE_(check_f<-24>(dilogarithm_f<1, 0, 2>(U_aphex{0.22, 0.22}, U_aphex{0.11, 0.11}), U_aphex{-0.048411922454435596,-0.04838741284417837}));
		TRUE_(check_f<-24>(dilogarithm_f<1, 0, 1>(U_aphex{0.22, 0.22}, U_aphex{0.11, 0.11}), U_aphex{-0.048411922454435596,-0.04838741284417837}));
		TRUE_(check_f<-36>(dilogarithm_f<1, 0, 0>(U_aphex{0.22, 0.22}, U_aphex{0.11, 0.11}), U_aphex{-0.048411922454435596,-0.04838741284417837}));

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
