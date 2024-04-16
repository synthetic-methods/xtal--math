#pragma once
#include "../any.hh"
#include "../dilute.hh"
#include "../discard.hh"
#include "../root.hh"
#include "../horner/all.hh"
#include "../pade/all.hh"

XTAL_ENV_(push)
namespace xtal::math::jacobi
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ord=0> struct zeta {};


////////////////////////////////////////////////////////////////////////////////

template <>
struct zeta<0>
{
	template <class S>
	class subtype: public S
	{
		using S_ = S;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_FN2 function(simplex_field_q auto &&t, complex_field_q auto &&s)
		XTAL_0EX
		{
			return function<N_lim>(_std::complex{XTAL_REF_(t)}, XTAL_REF_(s));
		}
		template <int N_lim=-1>
		XTAL_FN2 function(complex_field_q auto &&t, complex_field_q auto &&s)
		XTAL_0EX
		{
			unsigned constexpr M_lim = N_lim&0x7;
			using U_t = XTAL_TYP_(t); using V_t = typename U_t::value_type;
			using U_s = XTAL_TYP_(s); using V_s = typename U_s::value_type;
			using W_s = algebra::scalar_t<U_s[2]>;
			using M_s = algebra::series_t<W_s[2 + 2*M_lim]>;

			using re = bond::realize<U_s>;

			auto const e = _std::exp(s.imag() * -re::patio_1);
			auto const o = pade::unity_t<1, dilute<1>>::template function<4>(s.real());

			return function<N_lim>(XTAL_REF_(t), XTAL_REF_(s), M_s(o, e));
		}
		template <int N_lim=-1>
		XTAL_FN2 function(complex_field_q auto const &t, complex_field_q auto const &s, algebra::series_q auto const &oe)
		XTAL_0EX
		{
			unsigned constexpr M_lim = N_lim&0x7;
			using U_q = typename XTAL_TYP_(oe)::value_type::value_type;
			using U_s = XTAL_TYP_(s); using V_s = typename U_s::value_type;
			using U_t = XTAL_TYP_(t); using V_t = typename U_t::value_type;
		//	static_assert(is_q<U_t, U_s, U_q>);

			auto constexpr U_0 = U_s{ };
			auto constexpr V_1 = V_s{1};

			namespace diff_ = algebra::differential;
			using re = bond::realize<U_t>;
			using U_alpha = typename re::alpha_t;
			using U_sigma = typename re::sigma_t;
			using U_delta = typename re::delta_t;

			if constexpr (diff_::modular_q<V_t>) {
				auto  s0_im = -s.imag();
				auto  s0_re = -s.real();
				auto  const &[t2_re, t2_im] = reinterpret_cast<const V_t(&)[2]>(t);
				auto &t1_im = reinterpret_cast<const diff_::modular_t<iteratee_t<V_t>[1]> &>(t2_im);
				auto &t1_re = reinterpret_cast<const diff_::modular_t<iteratee_t<V_t>[1]> &>(t2_re);
			//	auto  t0_im = t1_im(0);
				auto  t0_re = t1_re(0);
				auto  w0_im = (t0_re*s0_im);
				auto  w0_re = (t1_re*s0_re + t1_im)(0);

				auto [x, y] = pade::wnity_t<1>::template function<4>(w0_re, w0_im).template reflected(1);
				
				return [&, x = XTAL_MOV_(x)] <size_t ...I>(bond::seek_t<I...>)
					XTAL_0FN_(U_0 +...+ (V_1/oe.get(1 + 2*I).sum(x)))
				(bond::seek_s<M_lim>{})*y - t0_re;
			}
			else {
				auto const &[s0_re, s0_im] = reinterpret_cast<const V_s(&)[2]>(s);
				auto  t0_re = t.real(); t0_re -= _std::round(t0_re);
				auto  t0_im = t.imag(); t0_im -= _std::round(t0_im);
				auto [x, y] = pade::wnity_t<1>::template function<4>(horner::term_f<-1>(t0_im, t0_re, s)).template reflected(1);
				
				return [&, x = XTAL_MOV_(x)] <size_t ...I>(bond::seek_t<I...>)
					XTAL_0FN_(U_0 +...+ (V_1/oe.get(1 + 2*I).sum(x)))
				(bond::seek_s<M_lim>{})*y - t0_re;
			}
		}

	};
};
template <int M_ord=0> using  zeta_t = process::confined_t<zeta<M_ord>>;
template <int M_ord=0>
XTAL_FN2 zeta_f(auto &&...oo)
XTAL_0EX
{
	return zeta_t<0>::function(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
