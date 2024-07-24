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
		using _op = bond::operating;
		using U = typename _op::alpha_type;
		//\
		using SVF = confined_t<filter<>>;
		using SVF = confined_t<prewarped<>, filter<>>;
		//\
		using Z = processor::monomer_t<prewarped<>, SVF>;
		using Z = processor::monomer_t<SVF>;

		SVF svf{};
		svf <<= occur::sample_t<>{44100};
		
		U constexpr omega = 2*2*3*3*5*5*7;
		U constexpr    up = 1;
		U constexpr    dn = 0;

		U _LP0{};
		U _LP1{};

		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(up, omega, _op::aphex_0)); TRUE_(_LP0 < _LP1); _LP0 = _LP1;
		
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;
		_LP1 = get<0>(svf(dn, omega, _op::aphex_0)); TRUE_(_LP0 > _LP1); _LP0 = _LP1;

		_std::array<U, 4> u_{1, 2, 0, 1};
		_std::array<U, 4> f_{3, 3, 3, 3};

		//\
		auto z = Z::braced_f(u_, f_);
		auto z = Z::braced_f(processor::let_f(u_), processor::let_f(f_));
//		TRUE_(true);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
