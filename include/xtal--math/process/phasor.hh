#pragma once
#include "./any.hh"

#include "../occur/indent.hh"
#include "../atom/phason.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ..._s>	struct  phasor;
template <typename ..._s>	using   phasor_t = confined_t<phasor<_s...>>;
template <typename ..._s>	concept phasor_q = bond::tag_in_p<phasor, _s...>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Manages a truncated fixed-point unit differential.

May be updated via the `indent` mechanism.
Provides evaluation/update via succession/replacement.

\todo    Attach a parameter to govern suboctave scaling?
*/
template <vector_q A, typename ...As>
struct phasor<A, As...>
{
	static auto constexpr N = _xtd::extent_v<A>;

	using U_phason = atom::math::phason_t<A>;
	using revalue_type = typename U_phason::revalue_type;
	using invalue_type = typename U_phason::invalue_type;
	using devalue_type = typename U_phason::devalue_type;
	using V_phason = atom::math::phason_t<invalue_type[N]>;

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
		using U_ = typename S_::head_type;

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS
		using S_::self;
		using S_::head;

	//	NOTE: Defined in-case `refine_head` is bypassed...
		XTAL_FN0_(go) (XTAL_DEF_(return,inline,explicit) operator U_,
			[] (auto &&that) XTAL_0FN_(to) (XTAL_REF_(that).head()))

		XTAL_DEF_(return,inline,set)
		bias()
		noexcept -> auto
		{
			return S_::template bias<revalue_type>();
		}

	public:// FLOW

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion> requires in_v<N_ion, -1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			if (o.head() == 0) {
				get<0>(head()) = zero;
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}

	public:// OPERATE
	//	TODO: Use `occur::resample` to manage downsampling via integer multiplication.

		/*!
		\brief   Evaluation by succession.
		\todo    Override constructors to apply fractional `bias`.
		*/		
		template <auto ...Is> requires (0 == sizeof...(Is))
		XTAL_DEF_(return,inline,let)
		method()
		noexcept -> decltype(auto)
		{
			using resample_type = occur::resample_t<>;
			if constexpr (complete_q<typename S_::template head_t<resample_type>>) {
				auto const rate = S_::template head<resample_type>().rate();
				auto &phi = ingress();
				auto  psi = phi.scaled(rate);
				return egress(psi.template apply<[] XTAL_1FN_(call) (bond::pack_f)>());
			}
			else {
				return egress(ingress());
			}
		}
	//	/*!
	//	\brief   Evaluation by (possibly indented) replacement then succession.
	//	*/		
	//	template <auto ...Is> requires (0 == sizeof...(Is))
	//	XTAL_DEF_(return,inline,let)
	//	method(fixed_shaped_q auto &&a)
	//	noexcept -> decltype(auto)
	//	{
	//		static_assert(fixed_shaped<decltype(a)>::extent() <= N);
	//		(void) S_::template flux<+1>(XTAL_REF_(a));
	//		return method();
	//	}
		/*!
		\returns The current differential after setting the stored frequency to `omega`.
		*/
		template <int N_root=1>
		XTAL_DEF_(return,let)
		method(revalue_type omega)
		noexcept -> auto
		{
			using resample_type = occur::resample_t<>;
			if constexpr (complete_q<typename S_::template head_t<resample_type>>) {
				omega *= S_::template head<resample_type>().period();
			}
			auto &u_phi = S_::head();
			u_phi[1] = U_phason::devalue_f(omega); return ++u_phi;
		}
		/*!
		\returns The current differential after scaling the incoming `phi` by `co`.
		*/
		template <int N_root=1>
		XTAL_DEF_(return,let)
		method(U_phason phi, revalue_type co)
		noexcept -> auto
		//	requires same_q<U_phason, typename S_::template head_t<ordinal_constant_t<1>>>
		{
			auto &u_phi = S_::head();
			phi *= co; u_phi[1] = phi[1]; return ++u_phi;
		}

	protected:
		/*!
		\brief   Evaluation by succession.
		*/		
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
		using U_ = typename S_::head_type;

	public:// ACCESS
		using S_::S_;
		using S_::self;
		using S_::head;

	//	NOTE:	Defined in-case `refine_head` is bypassed...
		XTAL_FN0_(go) (XTAL_DEF_(return,inline,implicit) operator U_,
			[] (auto &&o) XTAL_0FN_(to) (XTAL_REF_(o).head()))

	public:// OPERATE
		/*!
		\returns The current differential after scaling the incoming `phi` by `co`.
		*/
		template <int N_root=1>
		XTAL_DEF_(return,let)
		method(U_phason phi, revalue_type co)
		noexcept -> auto
		requires same_q<U_phason, typename S_::template head_t<ordinal_constant_t<1>>>
		{
			static_assert(2 == U_phason::size());
			auto &u_phi = S_::template head<constant_t<1>>();
			auto &v_phi = S_::template head<constant_t<0>>();
		//	Calculates the deviation of `phi[0]` w.r.t. phi[1],
		//	using the difference in `phi[1]` to determine the threshold for reset.

			u_phi[1]  = phi[1]; ++u_phi; auto i_phi = bond::math::bit_sign_f<devalue_type>(u_phi[0] == phi[0]);
			u_phi[0]  = phi[0];     phi *= co;
			v_phi[1]  = phi[1]; ++v_phi;
			v_phi[0] &= ~i_phi;
			v_phi[0] |=  i_phi&u_phi[0];

			return v_phi;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct context<process::math::phasor<_s...>>
{
	using superkind = context<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
	
	public:
		using S_::S_;

		using scale_type = occur::inferred_t<union scale, float>;

		template <extent_type N_mask=1>
		struct dispatch
		{
			using superkind = provision::voiced<void
			,	typename T_::index_type::template dispatch<N_mask>
			>;
			template <class R>
			using subtype = bond::compose_s<R
			,	typename S_::template dispatch<N_mask>
			,	superkind
			>;

		};

	};
};
template <scalar_q A>
struct context<process::math::phasor<A>> : context<process::math::phasor<A[2]>>
{
};
template <>
struct context<process::math::phasor< >> : context<process::math::phasor<typename bond::fit<>::alpha_type>>
{
};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
