#pragma once
#include "./any.hh"

#include "./power.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_exp=1, int M_cut=0>
struct   root;

template <int M_exp=1, int M_cut=0>
using    root_t = process::confined_t<root<M_exp, M_cut>>;

template <int M_exp=1, int M_cut=0, auto N_lim=0b11>
XTAL_DEF_(return,inline,let)
root_f(auto &&z)
noexcept -> decltype(auto)
{
	static_assert(M_exp != 0);
	XTAL_IF0
//	XTAL_0IF_(consteval) {return root_t<M_exp, 0    >::template static_method<   ~0>(XTAL_REF_(z));}
	XTAL_0IF (M_cut <=0) {return root_t<M_exp, 0    >::template static_method<N_lim>(XTAL_REF_(z));}
	XTAL_0IF (0 < M_cut) {return root_t<M_exp, M_cut>::template static_method<N_lim>(XTAL_REF_(z));}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_exp, int M_cut>
struct root
{
	static int constexpr M_exp_sgn = sign_n<M_exp>;
	static int constexpr M_exp_mag = M_exp*M_exp_sgn;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0b11>
		XTAL_DEF_(return,inline,set)
		static_method(auto &&z)
		noexcept -> auto
		{
			using _fix = bond::fixture<decltype(z)>;
			auto constexpr I_lim = below_m<(1<<4), (unsigned) N_lim>;
			XTAL_IF0
			XTAL_0IF (integral_variable_q<decltype(z)>) {
				return static_method<I_lim>(_fix::alpha_f(XTAL_REF_(z)));
			}
			XTAL_0IF_(return) (evaluate<I_lim>(XTAL_REF_(z)))
			XTAL_0IF_(return) (evaluate<I_lim>(XTAL_REF_(z), constant_n<2>))
			XTAL_0IF_(return) (evaluate<I_lim>(XTAL_REF_(z), constant_n<3>))
			XTAL_0IF_(return) (evaluate<I_lim>(XTAL_REF_(z), constant_n<5>))
			XTAL_0IF_(return) (evaluate<I_lim>(XTAL_REF_(z), constant_n<7>))
			XTAL_0IF_(else) {return pow(XTAL_REF_(z), _fix::alpha_1/M_exp);}
		}

	protected:
		template <int I_lim>
		XTAL_DEF_(return,inline,set)
		evaluate(auto &&z, constant_q auto i_exp)
		noexcept -> XTAL_ALL_(z)
		requires (M_exp%i_exp == 0 and 1 != M_exp/i_exp)
		{
			return root_t<M_exp/-i_exp>::template static_method<I_lim>(root_t<-i_exp>::template static_method<I_lim>(XTAL_REF_(z)));
		}
		template <int I_lim> requires in_n<M_exp_mag, 1>
		XTAL_DEF_(return,inline,set)
		evaluate(auto &&z)
		noexcept -> XTAL_ALL_(z)
		{
			using _fix = bond::fixture<decltype(z)>;

			XTAL_IF0
			XTAL_0IF (M_exp ==  1) {
				return XTAL_REF_(z);
			}
			XTAL_0IF (M_cut <= 0) {
				return one/(XTAL_REF_(z));
			}
			XTAL_0IF (0 <  M_cut) {
				return one/(XTAL_REF_(z) + _fix::minilon_f(M_cut));
			}
		}
		template <int I_lim> requires in_n<M_exp_mag, 2>
		XTAL_DEF_(return,inline,set)
		evaluate(complex_variable_q auto z)
		noexcept -> XTAL_ALL_(z)
		{
			using _fix = bond::fixture<decltype(z)>;

			XTAL_IF0
			XTAL_0IF_(consteval) {
				z *= _fix::haplo_f(M_exp >> 1);
				auto const x_re = z.real();
				auto const x_im = z.imag();
				auto const x_a2 = accumulator_f(x_re*x_re, x_im, x_im);
				auto const x_a1 = root_t<2>::template static_method<I_lim>(x_a2);

				auto y_re = x_a1 + x_re, v_re = y_re;
				auto y_im = x_a1 - x_re, v_im = y_im;
				if constexpr (M_exp < 0) {
					v_re *= x_a2;
					v_im *= x_a2;
				}
				y_re *= root_t<-M_exp_mag, 1>::template static_method<I_lim>(v_re);
				y_im *= root_t<-M_exp_mag, 1>::template static_method<I_lim>(v_im);

				auto const y_im_sgn = M_exp_sgn*_xtd::copysign(_fix::alpha_1, x_im);
				return {y_re, y_im*y_im_sgn};
			}
			XTAL_0IF_(else) {
				return root_f<M_exp_sgn, M_cut>(sqrt(z));
			}
		}
		template <int I_lim> requires in_n<M_exp_mag, 2, 3, 5, 7, 9>
		XTAL_DEF_(return,inline,set)
		evaluate(real_variable_q auto z)
		noexcept -> XTAL_ALL_(z)
		{
			auto constexpr z_one = XTAL_ALL_(z){1};
			auto const     z_sig = _xtd::copysign(z_one, z); z *= z_sig;
			XTAL_IF0
			XTAL_0IF_(consteval) {
				return z_sig*approximate<I_lim>(z);
			}
			XTAL_0IF (M_exp_mag != 2) {
				return z_sig*approximate<I_lim>(z);
			}
			XTAL_0IF (M_exp_mag == 2) {
				return z_sig*root_f<M_exp_sgn, M_cut>(sqrt(z));
			}
		}


		template <int I_lim>
		XTAL_DEF_(return,inline,set)
		approximate(real_variable_q auto z)
		noexcept -> XTAL_ALL_(z)
		{
			XTAL_IF0
			XTAL_0IF (0 < M_exp) {return exfunction<I_lim>(z);}
			XTAL_0IF (M_exp < 0) {return infunction<I_lim>(z);}
		}
		template <int I_lim>
		XTAL_DEF_(return,inline,set)
		exfunction(real_variable_q auto z)
		noexcept -> XTAL_ALL_(z)
		{
			return z*power_f<M_exp_mag - 1>(infunction<I_lim>(z));
		}
		template <int I_lim>
		XTAL_DEF_(return,inline,set)
		infunction(real_variable_q auto z)
		noexcept -> XTAL_ALL_(z)
		{
			using _fix = bond::fixture<decltype(z)>;
			using Z_sigma = typename _fix::sigma_type;
			using Z_delta = typename _fix::delta_type;
			using Z_alpha = typename _fix::alpha_type;

			Z_alpha constexpr  m_1 = _fix::dnsilon_f(_fix::exponent.depth + M_exp_mag);// Experimental error term...
			Z_delta constexpr  M_1 = _xtd::bit_cast<Z_delta>(m_1);

			Z_delta constexpr  N = -M_exp_mag;
			Z_alpha constexpr  n = -M_exp_mag;
			Z_alpha constexpr _n =  one/n;
			Z_delta constexpr _N =  M_1/N;

			XTAL_IF0
			XTAL_0IF_(consteval) {
				z += M_exp_mag*_fix::minilon_f(M_cut);
			}
			Z_delta constexpr  K_    = (N - one)*_N;
			Z_alpha constexpr  k_    = (n - one)*_n;
			Z_alpha const      z_    =         z*_n;
			Z_alpha constexpr  h     =         half;

			auto y = _xtd::bit_cast<Z_alpha>(K_ + _xtd::bit_cast<Z_delta>(z)/N);
			
			XTAL_IF0
			XTAL_0IF_(consteval) {
				auto v = z;
				for (unsigned i{}; i < 0x10 and v != y; ++i) {
					y *= accumulator_f(k_, z_, power_f<M_exp_mag>(v = y));
				}
				{
					y /= accumulator_f(h, h, z*power_f<M_exp_mag>(v = y));
				}
			}
			XTAL_0IF_(else) {
				#pragma unroll
				for (unsigned i{}; i < I_lim; ++i) {
					y *= accumulator_f(k_, z_, power_f<M_exp_mag>(y));
				}
			}
			return y;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
