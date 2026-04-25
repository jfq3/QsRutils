# QsRutils 0.2.0

* `avg_alpha()` now uses a compiled C++ rarefaction routine (`rrarefy_cpp`) in place of `vegan::rrarefy()`, yielding roughly a 34% end-to-end speedup.
