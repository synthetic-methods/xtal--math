#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
Wrapper used to tunnel an existing type using `std::tuple`-based traversal.

\see `process::cross`.
*/
////////////////////////////////////////////////////////////////////////////////

template <typename     ..._s> struct   indent;
template <typename     ..._s> concept  indent_q = bond::tag_in_p<indent, _s...>;
template <class S, int ...Ns> using    indent_s = bond::compose_s<S, indent<ordinal_constant_t<Ns>...>>;

template <constant_q ...Ns>
struct indent<Ns...>
{
	template <class S> using indicated_t = bond::pack_item_t<S, Ns{}...>;
	template <class S> using indicated_s = bond::pack_item_s<S, Ns{}...>;
	
	using superkind = bond::compose<void
	,	confined<bond::tag<indent>, confer<Ns>...>
	,	bond::compose_t<conferred_t>
	,	bond::compose_t<indicated_t>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(bond::pack_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using W_ =     indicated_s<S >;
		using U_ =   initializer_t<W_>;

	protected:
		/*!
		\brief   Converts the argument to an element Constructs a scalar fragment for `u`,
		         handling element conversion if supported by the container.

		\todo    Use strong-`value_type`s to map between fractional and floating-point values?
		*/
		XTAL_DEF_(return,inline,let)
		inject_f(auto &&u)
		noexcept -> auto
		{
			if constexpr      (requires{W_::devalue_f(XTAL_REF_(u));}
				and different_q<decltype(W_::devalue_f), _std::identity>
			)	{
				return W_::devalue_f(XTAL_REF_(u));
			}
			else {
				return XTAL_REF_(u);
			}
		}

	public:
		using S_::S_;//NOTE: Inherited and respecialized!

		/*!
		\brief   Constructs a scalar fragment for `u`,
		         handling element conversion if supported by the container.

		\todo    Use strong-`value_type`s to map between fractional and floating-point values?
		*/
		XTAL_NEW_(explicit)
		subtype(U_ u)
		noexcept requires infungible_q<S_, U_>
		:	S_{inject_f(XTAL_MOV_(u))}
		{}
		XTAL_NEW_(implicit)
		subtype(_std::initializer_list<U_> u_)
		noexcept requires make_p<S_, _std::initializer_list<U_>>
		:	S_{u_}
		{}

		template <extent_type N_mask=1>
		struct incept
		{
			using superkind = bond::compose<flow::mask<N_mask>, defer<indicated_t<S>>>;

			template <class R> requires un_n<sizeof...(Ns)>
			class subtype : public bond::compose_s<R, superkind>
			{
				static_assert(flow::any_q<R>);
				using R_ = bond::compose_s<R, superkind>;
				
			public:
				using R_::R_;
				using R_::self;
				using R_::head;

				/*!
				\brief   Forwards the message upstream.
				*/
				template <signed N_ion>
				XTAL_DEF_(return,inline,let)
				fuse(auto &&o)
				noexcept -> signed
				{
					return R_::template fuse<N_ion>(XTAL_REF_(o));
				}
				/*!
				\brief   Updates the internal state at the given path.
				*/
				template <signed N_ion> requires in_n<N_ion, +1>
				XTAL_DEF_(return,let)
				fuse(indent_q auto &&o)
				noexcept -> signed
				{
					/*/
					auto &m = bond::pack_item_f(o.seek(), head());
					XTAL_ALL_(m) x(o);
					_std::swap(m, x);
					return m == x;
					/*/
					auto &m = bond::pack_item_f(o.seek(), head());
					using M = XTAL_ALL_(m);
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
