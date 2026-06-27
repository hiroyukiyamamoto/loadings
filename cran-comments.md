## Test environments

* local Windows 11, R 4.4.0
* GitHub Actions:
  * windows-latest, R release
  * macos-latest, R release
  * ubuntu-latest, R release
  * ubuntu-latest, R devel

## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps ... NOTE
  unable to verify current time

## Changes

* Removed the package dependency on `geigen`.
* Updated generalized eigenvalue calculations and one-sided KPCA approximation.
