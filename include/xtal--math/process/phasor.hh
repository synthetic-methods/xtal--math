#pragma once
#include "./any.hh"

#include "../atom/phason.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ..._s>	struct  phasor;
template <typename ..._s>	using   phasor_t = confined_t<phasor<_s...>>;
template <typename ..._s>	concept phasor_q = bond::tagged_with_p<phasor, _s...>;


////////////////////////////////////////////////////////////////////////////////
///\
Manages a truncated fixed-point unit differential like `phasor`, \
providing evaluation/update via succession/replacement. \

///\todo\
Accommodate `std::complex` `value_type`s? \

template <vector_q A, typename ...As>
struct phasor<A, As...>
{
	static auto constexpr N = _xtd::extent_v<A>;

	using U_phason = atom::math::phason_t<A>;
	using coordinate_type = typename U_phason::coordinate_type;
	using inordinate_type = typename U_phason::inordinate_type;
	using   ordinate_type = typename U_phason::  ordinate_type;
	using V_phason = atom::math::phason_t<inordinate_type[N]>;

	
	using semikind = bond::compose<void
	//\
	,	refer<U_phason>
	,	cell::_detail::refer_multiplicative_group<U_phason>
	,	typename occur::math::indent_s<U_phason>::template incept<>
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
		Use `occur::resample` to manage downsampling \
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
			using resample_type = occur::resample_t<>;
			///\todo\
			Override constructors to apply fractional `bias`. \
			
			if constexpr (S_::template head_q<resample_type>) {
				auto const rate = S_::template head<resample_type>().rate();
				auto &phi = ingress();
				auto  psi = phi.scaled(rate);
				return egress(psi.template apply<[] XTAL_1FN_(call) (bond::pack_f)>());
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

			auto &u_phi = S_::template head<1>();
			auto &u_psi = S_::template head<0>();

		//	Calculates the deviation of `phi[0]` w.r.t. phi[1], \
		//	using the difference in `phi[1]` to determine the threshold for reset. \

			auto  u_delta = u_phi; u_delta[1] -= phi[1];
			//\
			auto &v_delta = u_delta.template self<inordinate_type>(constant_t<N>{});
			auto &v_delta = reinterpret_cast<V_phason const &>(u_delta);
			auto  n_delta = bond::math::bit_ceiling_f(aspect_f<unsigned>(v_delta[1]));
			auto  i_delta = condition_f<ordinate_type>(v_delta[0] >> n_delta);
			
			u_phi = XTAL_MOV_(phi);

			u_phi *= co;
			u_psi[1]  = u_phi[1];
			u_psi.operator++();
			u_psi[0] &=        ~i_delta;
			u_psi[0] |= u_phi[0]&i_delta;

			return u_psi;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
