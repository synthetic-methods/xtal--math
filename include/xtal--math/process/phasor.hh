#pragma once
#include "./any.hh"

#include "../atom/phason.hh"
#include "../atom/fourier/series.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ..._s> struct   phasor;
template <typename ..._s> using    phasor_t = confined_t<phasor<_s...>>;
template <typename ..._s> concept  phasor_q = bond::tag_p<phasor, _s...>;


////////////////////////////////////////////////////////////////////////////////
///\
Manages a truncated fixed-point unit differential like `phasor`, \
providing evaluation/update via succession/replacement. \

///\todo\
Accommodate `std::complex` `value_type`s? \

template <vector_q A, typename ...As>
struct phasor<A, As...>
{
	using _fit = bond::fit<A>;
	using _phi = atom::math::phason<A>;
	using coordinate_type = typename _phi::coordinate_type;
	using inordinate_type = typename _phi::inordinate_type;
	using   ordinate_type = typename _phi::  ordinate_type;

	static auto constexpr N  = _std::extent_v<A>;

	using U_lepton = atom::math::fourier::series_t<coordinate_type[N]>;
	using U_phason = atom::math::phason_t<coordinate_type[N]>;
	using V_phason = atom::math::phason_t<inordinate_type[N]>;
	
	using semikind = bond::compose<void
	//\
	,	refer<U_phason>
	,	cell::_detail::refer_multiplicative_group<U_phason>
	,	typename flow::indent_s<U_phason>::template afflux<>
	,	As...
	>;
	using superkind = bond::compose<bond::tag<phasor>
	,	semikind
	,	provision::biased<constant_t<1>>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS
		using S_::self;
		using S_::head;

		///\note\
		This is defined in-case `refine_head` is bypassed...

		///\todo\
		...find a cleaner way to define the conversion, perhaps via `refer`?

		XTAL_FX4_(to) (XTAL_DEF_(implicit operator) typename S_::head_type(), head())
		
		XTAL_DEF_(return,inline,set)
		bias()
		noexcept -> auto
		{
			return S_::template bias<coordinate_type>();
		}

	public:// OPERATE
		///\todo\
		Use `occur::sample` to manage downsampling \
		e.g. by integer multiplication followed by normalization. \

		///\
		Evaluation by (possibly indented) replacement then succession. \
		
		template <auto ...Is> requires (0 == sizeof...(Is))
		XTAL_DEF_(return,inline,let)
		method(fixed_shaped_q auto &&a)
		noexcept -> decltype(auto)
		{
			static_assert(fixed_shaped<decltype(a)>::extent() <= N);
			(void) S_::template flux<+1>(XTAL_REF_(a));
			return method();
		}

		///\
		Evaluation by uccession. \
		
		template <auto ...Is> requires (0 == sizeof...(Is))
		XTAL_DEF_(return,inline,let)
		method()
		noexcept -> decltype(auto)
		{
			///\todo\
			Override constructors to apply fractional `bias`. \
			
			if constexpr (requires {S_::sampling().rate();}) {
				auto const rate = S_::sampling().rate();
				auto &phi = ingress();
				XTAL_IF0
				XTAL_0IF (N == 1) {
					//\
					return egress(bond::pack_f(phi(0)));
					return egress(phi(0));
				}
				XTAL_0IF (N == 2) {
					return egress(bond::pack_f(phi(0), phi(1)*(rate)));
				}
				XTAL_0IF_(else) {
					return egress(phi.template apply<[] XTAL_1FN_(function) (bond::pack_f)>()*U_lepton(rate));
				}
			}
			else {
				return egress(ingress());
			}
		}

	protected:
		///\
		Evaluation by succession. \
		
		XTAL_DEF_(return,inline,let)
		ingress()
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (1 == bias()) {return ++head();}
			XTAL_0IF_(else)        {return   head();}
		};
		template <class Y>
		XTAL_DEF_(return,inline,let)
		egress(Y &&y)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (0 == bias()) {return (void) ++head(), XTAL_REF_(y);}
			XTAL_0IF_(else)        {return                  XTAL_REF_(y);}
		};
		
	};
	template <class S> requires phasor_q<bond::compose_s<S, semikind>>
	class subtype<S> : public bond::compose_s<S, semikind>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S, semikind>;

	public:// ACCESS
		using S_::S_;
		using S_::self;
		using S_::head;

		///\note\
		This is defined in-case `refine_head` is bypassed...

		XTAL_FX4_(to) (XTAL_DEF_(implicit operator) typename S_::head_type(), head())
		
	public:// REEVALUATION
		///\returns the current differential after scaling the incoming `phi` by `co`. \

		///\todo\
		Supply `precision` and/or `subdivision` `attach`ments? \

		template <int N_root=1>
		XTAL_DEF_(return,let)
		method(U_phason phi, coordinate_type co)
		noexcept -> auto
			requires same_q<U_phason, typename S_::template head_t<ordinal_constant_t<1>>>
		{
			static_assert(2 == U_phason::size());

			auto &_phi = S_::template head<1>();
			auto &_psi = S_::template head<0>();

		//	Calculates the deviation of `phi[0]` w.r.t. phi[1], \
		//	using the difference in `phi[1]` to determine the threshold for reset. \

			auto  u_delta = _phi - phi; u_delta[0] += phi[1];
			auto &v_delta = reinterpret_cast<V_phason const &>(u_delta);
			auto  n_delta = bond::math::bit_ceiling_f(magnum_f(v_delta[1]));
			auto  i_delta = condition_f<ordinate_type>(v_delta[0] >> n_delta);

			_phi = XTAL_MOV_(phi);

			_phi *= co;
			_psi[1]  = _phi[1];
			_psi.operator++();
			_psi[0] &=        ~i_delta;
			_psi[0] |= _phi[0]&i_delta;

			return _psi;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
