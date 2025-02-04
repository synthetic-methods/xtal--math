#pragma once
#include "./any.cc"
#include "./filter.hh"// testing...

#include "./prewarped.hh"



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
		using SVF = confined_t<filter<>>;
		using SVF = prewarped_t<ordinal_constant_t<1>, filter<>>;
		//\
		using Z = processor::monomer_t<prewarped<ordinal_constant_t<1>>, SVF>;
		using Z = processor::monomer_t<SVF>;

		SVF svf{};
		svf <<= typename occur::sample_t<>{44100};
		svf <<= typename filter<>::   limit_type{0};
		svf <<= typename filter<>::   order_type{2};
		svf <<= typename filter<>::  select_type{0};
		svf <<= typename filter<>::topology_type{0};
	
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
