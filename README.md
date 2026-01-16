XTAL MATH
=========

# About

XTAL--MATH is a DSP library that utilizes template composition as the means to build higher-order components from lower-level functional and stateful processes.

At the namespace level, the `math` namespaces specialize the `xtal` folder structure, with `xtal--math/*/` mapping to `xtal::*::math::`.
Within each `math` subspace, components may be further organized by an associated author or mathematician: currently Dirichlet, Fourier, Pade, Taylor, Zavalishin.


### Overview

| Folder                                                                                       | Description                                                                    |
|----------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
|[`atom/dirichlet/`                 ](include/xtal--math/atom/dirichlet/)                      | Structures for Dirichlet/Kronecker character/symbol generation.                |
|[`atom/fourier/`                   ](include/xtal--math/atom/fourier/)                        | Structures for linear and circular convolution and transformation.             |
|[`atom/pade/`                      ](include/xtal--math/atom/pade/)                           | Structure  for multiplicative complex sign and  mutually reciprocal magnitude. |
|[`process/pade/  `                 ](include/xtal--math/process/taylor/)                      | Approximations of the circular  and hyperbolic functions.                      |
|[`process/taylor/`                 ](include/xtal--math/process/taylor/)                      | Approximations of the logarithm and associated functions.                      |
|[`process/zavalishin/`             ](include/xtal--math/process/zavalishin/)                  | Dynamic processing for filtering, envelope generation, and differentiation.    |


### Underview

| File                                                                                         | Description                                                                    |
|----------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
|[`atom/phason.hh`                  ](include/xtal--math/atom/phason.hh?ts=3)                  | Compact fractional phase/frequency pair.                                       |
|[`atom/dirichlet/symbol.hh`        ](include/xtal--math/atom/dirichlet/symbol.hh?ts=3)        | Dirichlet/Kronecker characters/symbols.                                        |
|[`atom/fourier/serial.hh`          ](include/xtal--math/atom/fourier/serial.hh?ts=3)          | Linear convolution.                                                            |
|[`atom/fourier/series.hh`          ](include/xtal--math/atom/fourier/series.hh?ts=3)          | Circular convolution and transformation.                                       |
|[`atom/fourier/uniplex.hh`         ](include/xtal--math/atom/pade/uniplex.hh?ts=3)            | Compact multiplicative complex sign and mutually reciprocal magnitude.         |
|[`process/multiplex.hh`            ](include/xtal--math/process/multiplex.hh?ts=3)            | Modulation matrix.                                                             |
|[`process/pade/unity.hh`           ](include/xtal--math/process/pade/unity.hh?ts=3)           | Complex sinusoid with period `2π`, figuratively equivalent to `1^# &`.         |
|[`process/taylor/logarithm.hh`     ](include/xtal--math/process/taylor/logarithm.hh?ts=3)     | Approximations of the     logarithm  `Log`.                                    |
|[`process/taylor/monologarithm.hh` ](include/xtal--math/process/taylor/monologarithm.hh?ts=3) | Approximations of the monologarithm `-Log[1 - #] &`.                           |
|[`process/zavalishin/curve.hh`     ](include/xtal--math/process/zavalishin/curve.hh?ts=3)     | (Anti)saturation curves.                                                       |
|[`process/zavalishin/differ.hh`    ](include/xtal--math/process/zavalishin/differ.hh?ts=3)    | First-order differences with frequency/amplitude normalization.                |
|[`process/zavalishin/filter.hh`    ](include/xtal--math/process/zavalishin/filter.hh?ts=3)    | Arbitrary-order State Variable Filter (SVF).                                   |
|[`process/zavalishin/vactrol.hh`   ](include/xtal--math/process/zavalishin/vactrol.hh?ts=3)   | Envelope generation using `filter`.                                            |


(Note, as with [`xtal`](///github.com/synthetic-methods/xtal),
each `*.hh` is accompanied by a test-suite `*.cc` which is not included by `CMake`/`Conan`.)


# Usage

At the class level, components are parameterized both by the containing `struct` and the enclosed `method`s.
The intermediate template `subtype` provides the mechanism for decorator composition via `xtal::bond::compose`.

```c++
template <auto ...Ms>
struct unit
{
   template <class S>
   class subtype
   {
      template <auto ...Ns>
      auto method(auto ...xs) -> noexcept decltype(auto)
      {
      // ...   
      }
   };
};
```

The outer `Ms...` of the `struct` are used to specify compile-time configuration.
The inner `Ns...` of the `method` are enumerated and made available via _vtable_.

```c++
template <int M_ism=1>
struct sine
{
   template <class S>
   class subtype
   {
      template <int N_ord=1>
      auto method(auto ...xs) -> noexcept decltype(auto)
      {
      // Evaluate polynomial of order `N_ord`...
      // ...when `M_ism ==  2` approximate the hyperbolic function.
      // ...when `M_ism ==  1` approximate the   circular function.
      // ...when `M_ism == -1` approximate the   circular arc.
      // ...when `M_ism == -2` approximate the hyperbolic area.
      }
   };
};
```

Instances of the `struct` may thus be composed with other components.

```c++
using haversine = bond::compose<square<1>, sine<1>, halve<1>>;
```

For more information on the composition mechanism, see the [`xtal`](///github.com/synthetic-methods/xtal) documentation.
Detail on individual components can be found in the comments or associated Doxygen output.

The most comprehensive reference for usage is currently the test suite (sorry!).
The following annotated example comes from the
[`zavalishin/filter.hh`](include/xtal--math/process/zavalishin/filter.hh?ts=3) test-suite.


# Example

The following test-code and output demonstrates a polyphonic ringing filter
configured as a polyphonic instrument.

## Code

```c++
// Process definition...

using T_content = process::filter<U_alpha[2]>;         // 2nd-order filter.
using T_context = occur::context_t<T_content>;         // 2nd-order filter parameters.

using P_damp    = typename T_context::damp_parameter;
using P_fade    = typename T_context::fade_parameter;
using Q_order   = typename T_context::order_attribute;
using U_stage   = typename T_context::stage_type;      // Note stage: `0` on, `1` off, `-1` cut/rest.
using U_event   = flow::key_s<U_stage>;                // Key-Trigger pair.

// Process composition...

using T_process = process::confined_t<void             // Wrap the process.
,  provision::math::prewarping<0>                      // Apply s-plane frequency-prewarping to the first argument.
,  pulse< 1>                                           // Fulfil first argument with gate-signal controlled by `stage`.
,  reuse< 0>                                           // Resets the filter-state when note is initialized.
,  reuse<-1>                                           // Provides release detection.
,  typename U_stage   ::template   assignment<P_damp>  // Creates a table associating `stage` with damping.
,  typename P_damp    ::template   attend<>            // Attach/append damping to the arguments-list.
,  typename P_fade    ::template   attend<>            // Attach/append  fading to the arguments-list.
,  typename T_context ::template   attach<>            // Attach the remaining   object properties.
,  typename T_context ::template dispatch<>            // Attach the remaining template parameters.
,  T_content
>;

// Processor composition...

using T_scheduler = schedule::slicer_t<provision::spooled<extent_constant_t<0x10>>>;

using T_processor = processor::polymer_t<T_process     // Map the process while applying with polyphony.
,  T_scheduler::template accept<U_event>               // Moderate control via buffer-slicing scheduler.
,  provision::stored <null_type[0x100]>                // Use `std::vector`-based buffering.
,  provision::spooled<null_type[0x100]>                // Use `std::vector`-based event-spooliing.
>;

auto z_resize =   occur::resize_t<>(0x020);            // Let buffer-length equal `16`.
auto z_cursor =   occur::cursor_t<>(0x020);            // Let window-length equal `16`.
auto z_sample =   occur::resample_f(2*2*3*3*5*5*7*7);  // Let sampling-rate    equal `44100 Hz`.
auto z_omega  =   processor::let_f(59*61);             // Let filter-frequency equal `3599  Hz`.
auto z        = T_processor::bind_f(z_omega);

z <<= Q_order {2};                                     // Use the second-order filter (template).
z <<= P_damp{1};                                       // Initialize with maximum-damping.
z <<= P_fade{1};                                       // Initialize with maximum-volume.
z <<= flow::assign_f(U_stage{ 0}) << P_damp{0.000F};   // Associate note-on  with maximum-resonance.
z <<= flow::assign_f(U_stage{ 1}) << P_damp{0.060F};   // Associate note-off with  medium-resonance.
z <<= flow::assign_f(U_stage{-1}) << P_damp{0.707F};   // Associate note-cut with minimal-resonance.

z <<= z_sample;                                        // Influx the current sample-rate.
z <<= z_resize;                                        // Influx the current buffer-size.

z.lead() >>= typename T_context::stage_type{-1};       // Efflux all voices to note-cut.

z >>= flow::cue_f(0x08).then(U_event{69, 0});          // Cue note-on  for `A4` at sample `0x08`.
z >>= flow::cue_f(0x18).then(U_event{69, 0});          // Cue note-on  for `A4` at sample `0x18`, cutting the previous note.
z >>= flow::cue_f(0x28).then(U_event{69, 1});          // Cue note-off for `A4` at sample `0x08` within next block.

assert(0 == z.ensemble().size());                      // Expect `0` voices currently allocated.
assert(0 == z.efflux(z_cursor++));                     // Expect successful render!
plot<27>(z.store(), 0x08, 0x18);

z >>= flow::cue_f(0x20).then(U_event{69, 1});          // Cue note-cut for `A4` at sample `0x20`.

assert(2 >= z.ensemble().size());                      // Expect `2` voices currently decaying.
assert(0 == z.efflux(z_cursor++));                     // Expect successful render!
plot<27>(z.store(), 0x08);

assert(1 >= z.ensemble().size());                      // Expect `1` voices currently decaying.
assert(0 == z.efflux(z_cursor++));                     // Expect successful render!
plot<27>(z.store(), 0x00);
```

## Plot

The test-output is presented below, annotated with the note events marked on the central axis.

```
────────────────────────────────────────────────────────
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                           ⊣━━━━━━                       A4 note-on
                            ━━━━━━━━━━━━━━━━━━           
                            ━━━━━━━━━━━━━━━━━━━━━━━━━    
                            ━━━━━━━━━━━━━━━━━━━━━━━━━╸   
                            ━━━━━━━━━━━━━━━━━━━╸         
                            ━━━━━━━━                     
                       ╺━━━━                             
           ━━━━━━━━━━━━━━━━━                             
   ╺━━━━━━━━━━━━━━━━━━━━━━━━                             
  ━━━━━━━━━━━━━━━━━━━━━━━━━━                             
       ╺━━━━━━━━━━━━━━━━━━━━                             
                  ╺━━━━━━━━━                             
                            ━━╸                          
                            ━━━━━━━━━━━━━━━╸             
                            ━━━━━━━━━━━━━━━━━━━━━━━━     
                            ━━━━━━━━━━━━━━━━━━━━━━━━━━   
                           ⊣━━━━━━━━━━━━━━━━━╸           A4 note-on/cut
                            ━━━━━━━╸                     
                            ━━━━━━╸                      
                            ━━━━━━                       
                            ━━━╸                         
                          ━━                             
                 ━━━━━━━━━━━                             
        ━━━━━━━━━━━━━━━━━━━━                             
┄┄╴╺━━━━━━━━━━━━━━━━━━━━━━━━ ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄ 
  ━━━━━━━━━━━━━━━━━━━━━━━━━━                             
       ╺━━━━━━━━━━━━━━━━━━━━                             
                  ╺━━━━━━━━━                             
                            ━━╸                          
                            ━━━━━━━━━━━━━━━╸             
                            ━━━━━━━━━━━━━━━━━━━━━━━━     
                            ━━━━━━━━━━━━━━━━━━━━━━━━━━   
                           ⊣━━━━━━━━━━━━━━               A4 note-off
                     ━━━━━━━                             
  ━━━━━━━━━━━━━━━━━━━━━━━━━━                             
━━━━━━━━━━━━━━━━━━━━━━━━━━━━                             
━━━━━━━━━━━━━━━━━━━━━━━━━━━━                             
━━━━━━━━━━━━━━━━━━━━━━━━━━━━                             
              ━━━━━━━━━━━━━━                             
                            ━━━                          
                            ━━━━━━━━━━━━━━━━━━━━         
                            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 
                            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 
                            ━━━━━━━━━━━━━━━━━━━━━━━━━━╸  
                            ━━━━━━━━━━━━━━               
                           ╺                             
             ━━━━━━━━━━━━━━━                             
   ╺━━━━━━━━━━━━━━━━━━━━━━━━                             
╺━━━━━━━━━━━━━━━━━━━━━━━━━━━                             
     ━━━━━━━━━━━━━━━━━━━━━━━                             
               ━━━━━━━━━━━━━                             
                                                         
                            ━━━━━━━━━━╸                  
                            ━━━━━━━━━━━━━━━━━━━╸         
                            ━━━━━━━━━━━━━━━━━━━━━━━      
                            ━━━━━━━━━━━━━━━━━━━━         
┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄⊣━━━━━━━━━ ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄ A4 note-cut
                          ╺━                             
                     ━━━━━━━                             
                    ━━━━━━━━                             
                     ━━━━━━━                             
                       ━━━━━                             
                         ━━━                             
                           ━                             
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
────────────────────────────────────────────────────────
