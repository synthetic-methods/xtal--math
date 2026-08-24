#pragma once
#include "./any.hh"

#include <initializer_list>




XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
Wrapper used to tunnel an existing type using `std::tuple`-based traversal.

\see `process::cross`.
*/
////////////////////////////////////////////////////////////////////////////////

template <             typename ..._s> XTAL_TYP_(new)    dent;
template <class T                    > XTAL_TYP_(new) in_dent    {using type =          T           ;};
template <bond::tag_inner_q<dent>  T > XTAL_TYP_(new) in_dent<T> {using type = typename T::data_type;};
template <class T                    > XTAL_TYP_(set) in_dent_t = typename in_dent<based_t<T>>::type;
template <             typename ...Ts> XTAL_TYP_(ask) in_dent_q = bond::tag_inner_p<dent, Ts...>;
template <             typename ...Ts> XTAL_TYP_(ask)    dent_q = in_v<true, in_dent_q<Ts>...> and xtd::fungible_with<in_dent_t<Ts>...>;
template <class S,     int      ...Ns> XTAL_TYP_(set)    dent_s = bond::compose_s<S, dent<ordinal_constant_t<Ns>...>>;

template <constant_q ...Ns>
struct dent<Ns...>
{
	template <class S> using conferred_t =        confined_t<defer<S>>;
	template <class S> using indicated_t = bond::pack_item_t<S, Ns{}...>;
	template <class S> using indicated_s = bond::pack_item_s<S, Ns{}...>;

	using superkind = bond::compose<void
	,	confined<bond::tag<dent>, confer<Ns>...>
	,	bond::compose_t<conferred_t>
	,	bond::compose_t<indicated_t>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(bond::pack_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using C_ =     indicated_s<S >;
		using X_ =     indicated_t<S >;

		XTAL_VAL_(set) use_devalue = requires {
			{C_::devalue_f} -> different_q<std::identity>;
		};

	public:
		/*/
		using S_::S_;
		/*/
		XTAL_VAL_(delete) (subtype, noexcept=default)
		XTAL_VAL_(create) (subtype, noexcept=default)
		XTAL_VAL_(move)   (subtype, noexcept=default)
		XTAL_VAL_(copy)   (subtype, noexcept=default)
		XTAL_VAL_(induce) (subtype, noexcept,subtype)
	//	XTAL_VAL_(reduce) (subtype, noexcept,S_)
		/***/
		using data_type = S;

		XTAL_VAL_(new,explicit)
		subtype(X_ x_)
		noexcept
		:	S_(XTAL_MOV_(x_))
		{}
		/*!
		\brief   Constructs a fragment, using conversion as supported by the container.
		\todo    Use strong-`value_type`s to map between fractional and floating-point values?
		*/
		XTAL_VAL_(new,explicit)
		subtype(auto &&...oo)
		noexcept
		requires un_v<use_devalue> and requires  {X_{              XTAL_REF_(oo) ...};} and
		requires {X_{              XTAL_REF_(oo) ...};}
		:	    S_(X_{              XTAL_REF_(oo) ...} )
		{}
		XTAL_VAL_(new,explicit)
		subtype(auto &&...oo)
		noexcept
		requires in_v<use_devalue> and requires  {X_{C_::devalue_f(XTAL_REF_(oo))...};} and
		requires {X_{C_::devalue_f(XTAL_REF_(oo))...};}
		:	    S_(X_{C_::devalue_f(XTAL_REF_(oo))...} )
		{}

		template <extent_type N_mask=1>
		struct attach
		{
			using superkind = bond::compose<flow::mask<N_mask>, defer<S>>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
			//	static_assert(0 == sizeof...(Ns));
				static_assert(flow::any_q<R>);
				using R_ = bond::compose_s<R, superkind>;

			public:
				using R_::R_;
				using R_::self;
				using R_::head;
				using R_::headed;

				/*!
				\brief   Forwards the message upstream.
				*/
				template <signed N_ion>
				XTAL_VAL_(return,inline,let)
				fuse(auto &&o)
				noexcept -> signed
				{
					return R_::template fuse<N_ion>(XTAL_REF_(o));
				}
				/*!
				\brief   Updates the internal state at the given path.
				*/
				template <signed N_ion>// requires in_v<N_ion, +1>
				XTAL_VAL_(return,inline,let)
				fuse(dent_q<S> auto &&o)
				noexcept -> signed
				{
					auto &m = bond::pack_item_f(o.seek(), head());
					using M = XTAL_ALL_(m);
					/*/
					M x(o); std::swap(m, x);
					return m == x;
					/*/
					m.~M(); new (&m) M(XTAL_REF_(o));
					return 0;
					/***/
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
