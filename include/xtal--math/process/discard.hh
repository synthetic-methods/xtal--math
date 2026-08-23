#pragma once
#include "./any.hh"

#include "./root.hh"
#include "./identity.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_car=0, int M_aux=0>
struct  discard;

template <int M_car=0, int M_aux=0>
using   discard_t = process::confined_t<discard<M_car, M_aux>>;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_car, int M_aux>
struct discard
{
//	TODO: Handle `M_aux` for the base-class.

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_car=0>
		XTAL_VAL_(return,inline,let)
		method(auto const &x) const
		noexcept -> auto
		{
			auto constexpr m_car = M_car&1;
			auto constexpr n_car = N_car&1;
			XTAL_IF0
			XTAL_0IF (m_car == n_car) {return XTAL_ALL_(x) {one};}
			XTAL_0IF (m_car ==     0) {return                 x ;}
			XTAL_0IF (m_car ==     1) {return   root_f<-1, 1>(x);}
		}
		template <int N_car=0>
		XTAL_VAL_(return,inline,let)
		method(auto const &x, auto const &z) const
		noexcept -> XTAL_ALL_(x)
		{
			auto constexpr m_car = M_car&1;
			auto constexpr n_car = N_car&1;
			XTAL_IF0
			XTAL_0IF (m_car == n_car) {return XTAL_ALL_(x) {one};}
			XTAL_0IF (m_car ==     0) {return                x  ;}
			XTAL_0IF (m_car ==     1) {return root_f<-1, 1>(x*z);}
		}

	};
};
template <int M_aux>
struct discard<0, M_aux>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int M_arg=0>
		struct transfix
		{
		//	using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			using subtype = bond::compose_s<R>;

		};

	};
};
template <int M_aux>
struct discard<1, M_aux>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int M_arg=0>
		struct transfix
		{
			static int constexpr M_pow = signum_v<M_aux, 1>;

		//	using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:
				using R_::R_;

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo) const
				noexcept -> decltype(auto)
				requires
				requires (R_ const &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo)
				noexcept -> decltype(auto)
				requires
				requires (R_       &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}

			private:
				template <auto ...Ns> requires (0 == M_arg)
				XTAL_VAL_(return,inline,set)
				method_f(auto &&_, auto &&u, auto &&...oo)
				noexcept -> decltype(auto)
				{
					auto &r_ = xtd::qualify_cast<R_>(XTAL_REF_(_));
					auto  v  = r_.template method<Ns...>(u, XTAL_REF_(oo)...);
					using V  = XTAL_ALL_(v);
					using U  = XTAL_ALL_(u);
					static_assert(same_q<U, V>);
					return v*root_f<M_pow, 1>(XTAL_REF_(u));
				}
				template <auto ...Ns> requires (0 != M_arg)
				XTAL_VAL_(return,inline,set)
				method_f(auto &&_,           auto &&...ou)
				noexcept -> decltype(auto)
				{
					auto &r_ = xtd::qualify_cast<R_>(XTAL_REF_(_));
					auto  u  = bond::seek_argument_f<M_arg>(ou...);
					auto  v  = r_.template method<Ns...>(   XTAL_MOV_(ou)...);
					using U  = XTAL_ALL_(u);
					using V  = XTAL_ALL_(v);
					static_assert(same_q<U, V>);
					return v*root_f<M_pow, 1>(XTAL_MOV_(u));
				}

			};
		};

	};
};
template <>
struct discard<1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int M_arg=0>
		struct transfix
		{
		//	using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:
				using R_::R_;

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo) const
				noexcept -> decltype(auto)
				requires
				requires (R_ const &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo)
				noexcept -> decltype(auto)
				requires
				requires (R_       &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}
			
			private:
				template <auto ...Ns> requires (0 == M_arg)
				XTAL_VAL_(return,inline,set)
				method_f(auto &&_, auto &&u, auto &&...oo)
				noexcept -> auto
				{
					auto &r_ = xtd::qualify_cast<R_>(XTAL_REF_(_));
					auto  v  = r_.template method<Ns...>(u, XTAL_REF_(oo)...);
					using V  = XTAL_ALL_(v);
					using U  = XTAL_ALL_(u);
					XTAL_IF0
					XTAL_0IF (same_q<U, V>) {
						return v*XTAL_REF_(u);
					}
					XTAL_0IF (complex_variable_q<V>) {
						XTAL_IF1_(consteval) {
							return V{v.real(), v.imag()*XTAL_REF_(u)};
						}
						XTAL_0IF_(else) {
							destruct_f<1>(v) *= XTAL_REF_(u);
							return v;
						}
					}
					XTAL_0IF (complex_field_q<V>) {
						return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
					}
				}
				template <auto ...Ns> requires (0 != M_arg)
				XTAL_VAL_(return,inline,set)
				method_f(auto &&_, auto ...ou)
				noexcept -> auto
				{
					auto &r_ = xtd::qualify_cast<R_>(XTAL_REF_(_));
					auto  u  = bond::seek_argument_f<M_arg>(ou...);
					auto  v  = r_.template method<Ns...>(   XTAL_MOV_(ou)...);
					using U  = XTAL_ALL_(u);
					using V  = XTAL_ALL_(v);
					XTAL_IF0
					XTAL_0IF (same_q<U, V>) {
						return v*XTAL_MOV_(u);
					}
					XTAL_0IF (complex_variable_q<V>) {
						XTAL_IF1_(consteval) {
							return V{v.real(), v.imag()*XTAL_MOV_(u)};
						}
						XTAL_0IF_(else) {
							destruct_f<1>(v) *= XTAL_MOV_(u);
							return v;
						}
					}
					XTAL_0IF (complex_field_q<V>) {
						return complexion_f(v.real(), v.imag()*XTAL_MOV_(u));
					}
				}

			};
		};

	};
};
template <int M_aux>
struct discard<2, M_aux>
{
private:
	XTAL_VAL_(set) I_aux = sigpar_v<M_aux>;

public:
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int M_arg=0>
		struct transfix
		{
		//	using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:
				using R_::R_;

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo) const
				noexcept -> decltype(auto)
				requires
				requires (R_ const &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo)
				noexcept -> decltype(auto)
				requires
				requires (R_       &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}

			private:
				template <auto ...Ns> requires (0 == M_arg)
				XTAL_VAL_(return,inline,set)
				method_f(auto &&_, auto &&u, auto &&...oo)
				noexcept -> decltype(auto)
				{
					auto constexpr i = unstruct_t<decltype(u)>(I_aux);
					return xtd::qualify_cast<R_>(XTAL_REF_(_)).
						template method<Ns...>(i*square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
				}
				template <auto ...Ns> requires (0 != M_arg)
				XTAL_VAL_(return,inline,set)
				method_f(auto &&_,             auto ...ou)
				noexcept -> decltype(auto)
				{
					auto         & u = bond::seek_argument_f<M_arg>(ou...);
					auto constexpr i = unstruct_t<decltype(u)>(I_aux);
					u *= u;
					u *= i;
					return xtd::qualify_cast<R_>(XTAL_REF_(_)).
						template method<Ns...>(XTAL_MOV_(ou)...);
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
