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
		XTAL_DEF_(return,inline,set)
		method_f(auto const &x)
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
		XTAL_DEF_(return,inline,set)
		method_f(auto const &x, auto const &z)
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
		struct infix
		{
			using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			using subtype = bond::compose_s<R, superkind>;

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
		struct infix
		{
			static int constexpr M_pow = sign_v<M_aux, 1>;

			using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;

			public:
				using R_::R_;

			//	TODO: Account for `const &` and `&`, or find a cleaner way to express...
				XTAL_FN2_(do) (template <auto ...Is>
				XTAL_DEF_(return,inline,let)
				method(auto &&u, auto &&...oo),
				noexcept -> auto
				{
					auto  v = R_::template method<Is...>(u, XTAL_REF_(oo)...);
					using V = XTAL_ALL_(v);
					using U = XTAL_ALL_(u);
					static_assert(same_q<U, V>);

					return v*root_f<M_pow, 1>(XTAL_REF_(u));
				})
				template <auto ...Is>
				XTAL_DEF_(return,inline,set)
				method_f(auto &&u, auto &&...oo)
				noexcept -> auto
				{
					auto  v = R_::template method_f<Is...>(u, XTAL_REF_(oo)...);
					using V = XTAL_ALL_(v);
					using U = XTAL_ALL_(u);
					static_assert(same_q<U, V>);

					return v*root_f<M_pow, 1>(XTAL_REF_(u));
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
		struct infix
		{
			using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;

			public:
				using R_::R_;

				XTAL_FN2_(do) (template <auto ...Is>
				XTAL_DEF_(return,inline,let)
				method(auto &&u, auto &&...oo),
				noexcept -> auto
				{
					auto  v = R_::template method<Is...>(u, XTAL_REF_(oo)...);
					using V = XTAL_ALL_(v);
					using U = XTAL_ALL_(u);

					XTAL_IF0
					XTAL_0IF (same_q<U, V>) {
						return v*XTAL_REF_(u);
					}
					XTAL_0IF (complex_variable_q<V>) {
						XTAL_IF1_(consteval) {
							return V{v.real(), v.imag()*XTAL_REF_(u)};
						}
						XTAL_0IF_(else) {
							destruct_f(v)[1] *= XTAL_REF_(u); return v;
						}
					}
					XTAL_0IF (complex_field_q<V>) {
						return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
					}
				})
				template <auto ...Is>
				XTAL_DEF_(return,inline,set)
				method_f(auto &&u, auto &&...oo)
				noexcept -> auto
				{
					auto  v = R_::template method_f<Is...>(u, XTAL_REF_(oo)...);
					using V = XTAL_ALL_(v);
					using U = XTAL_ALL_(u);

					XTAL_IF0
					XTAL_0IF (same_q<U, V>) {
						return v*XTAL_REF_(u);
					}
					XTAL_0IF (complex_variable_q<V>) {
						XTAL_IF1_(consteval) {
							return V{v.real(), v.imag()*XTAL_REF_(u)};
						}
						XTAL_0IF_(else) {
							destruct_f(v)[1] *= XTAL_REF_(u); return v;
						}
					}
					XTAL_0IF (complex_field_q<V>) {
						return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
					}
				}

			};
		};

	};
};
template <int M_aux>
struct discard<2, M_aux>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int M_arg=0>
		struct infix
		{
			using superkind = typename identity_t<>::template infix<M_arg>;
			
			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;

			public:
				using R_::R_;

				XTAL_FN2_(do) (template <auto ...Is>
				XTAL_DEF_(return,inline,let)
				method(auto &&u, auto &&...oo),
				noexcept -> decltype(auto)
				{
					using _fit = bond::fit<decltype(u)>;
					auto constexpr v = _fit::alpha_f(cosign_v<M_aux>);
					return R_::template        method<Is...>(v*square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
				})
				template <auto ...Is>
				XTAL_DEF_(return,inline,set)
				method_f(auto &&u, auto &&...oo)
				noexcept -> decltype(auto)
				{
					using _fit = bond::fit<decltype(u)>;
					auto constexpr v = _fit::alpha_f(cosign_v<M_aux>);
					return R_::template method_f<Is...>(v*square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
