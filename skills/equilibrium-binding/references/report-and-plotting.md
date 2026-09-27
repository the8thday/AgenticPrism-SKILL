# Rendering and reporting

`agentic-prism render --run PATH --style prism_like` consumes immutable saved
results. It checks scientific artifact hashes, writes both themes and a separate
render manifest, and never invokes fitting. An altered input/result/config is an
error; create a new analysis run to change scientific settings.

Both themes share points, predictions, axis transforms/ranges and group order.
Style tokens live in the shared renderer. Prism-like uses white background,
left/bottom axes, outward ticks, clear markers and restrained color. It is not an
official GraphPad template. SVG and PDF are vectors; PNG is 300 or 600 dpi. Figure
widths are exactly 85 or 180 mm; height is 1.2 times width. Font fallback is logged.
SVG glyphs are paths and PDF embeds supported fonts to improve portability.

The zero-concentration control is in a separate small linear panel, never moved
to a fictitious positive log coordinate. All included observations are visible;
explicit exclusions are red crosses. Same-curve, same-concentration repeated
observations have descriptive mean ± SD bars when n > 1. They are not an
independent-experiment CI or a confidence band; fitted values still use the
original observations. Cross-curve/experiment error bars are not produced.

Offline HTML embeds all images and all download content (including PDF/PNG).
The theme selector changes visible images and their corresponding downloads.
No CDN, network fetch, refit or external display resources are required.
Public examples keep source attribution. A report includes failures and limited
curves, status, intervals, baseline/amplitude, methods, configuration, sensitivity,
independent-experiment summaries, artifact downloads and a rerun command.

Before delivery inspect actual rendered figures and HTML. Check long IDs/units,
fonts, legends, mobile width, theme switching, downloads and offline behavior.
Neither file existence nor structural Skill validation proves scientific accuracy.
