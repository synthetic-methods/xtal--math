#pragma once
#include "./any.cc"

#include "./prewarped.hh"
#include "./staged.hh"


#include "./filter.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("filter")
{
	TRY_("instantiation")
	{
		using _fit = bond::fit<>;
		using U_alpha = typename _fit::alpha_type;
		//\
		using SVF = confined_t<staged<0>, filter<>>;
		using SVF = prewarped_t<ordinal_constant_t<1>, staged<0>, filter<>>;

		using meta = any<filter<>>;
		//\
		using Z = processor::monomer_t<prewarped<ordinal_constant_t<1>>, SVF>;
		using Z = processor::monomer_t<SVF>;

		SVF svf{};
		svf <<= typename occur::sampling_t<>{44100};
		svf <<= typename meta::   limit_type{0};
		svf <<= typename meta::   order_type{2};
		svf <<= typename meta::  select_type{0};
		svf <<= typename meta::patch_type{0};
	
		U_alpha constexpr omega = 2*2*3*3*5*5*7;
		U_alpha constexpr   rho = 1;
		U_alpha constexpr    up = 1;
		U_alpha constexpr    dn = 0;

		U_alpha _LP0{};
		U_alpha _LP1{};

		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = svf(up, omega, rho); TRUE_(_LP0 < _LP1); _LP0 = _LP1;

		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = svf(dn, omega, rho); TRUE_(_LP0 > _LP1); _LP0 = _LP1;

		_std::array<U_alpha, 4> u_{1, 2, 0, 1};
		_std::array<U_alpha, 4> f_{3, 3, 3, 3};

		//\
		auto z = Z::bind_f(u_, f_);
//		auto z = Z::bind_f(processor::let_f(u_), processor::let_f(f_));
//		TRUE_(true);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
