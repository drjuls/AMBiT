## How to contribute to AMBiT
Thank you for considering contributing to AMBiT! AMBiT is under constant development and we welcome 
contributions and suggestions to improve the software. In particular, we encourage bug reports and 
bug-fixes, new features and improvements to existing features, and improvements to documentation.

Before you start contributing, though, take a moment to look over the following guidelines to make sure
the everything runs as smoothly as possible. 

### Reporting a bug
Rather than emailing the developers, please run through this checklist to file a bug report:

1. Make sure you're using the latest version of AMBiT (on either the `dev` or `master` branch) - the bug
   may have already been fixed in a more recent update.
2. Check the list of open [issues](https://github.com/drjuls/AMBiT/issues) to see if anyone has
   already reported this. If there's already an open issue on this bug, then please add a comment if 
   you have any new information to add (e.g. more information on how to trigger the bug, or other 
   circumstances it arises in).
3. If no-one has reported the bug yet, [open a new issue](https://github.com/drjuls/AMBiT/issues/new),
   with a clear description of the unintended behaviour and how it differs from what you expected. 
   Please provide as much relevant information in the description as possible. In particular, make sure 
   to include:
    
    * The input file which triggers the unexpected behaviour, as well as the resulting output. Ideally,
      try and remove everything from the input which is not necessary to trigger the bug (this is called
      a minimal working example). Make a [GitHub Gist](https://gist.github.com/) for these files,
      rather than pasting them directly into the issue.
    * The version of AMBiT you are running, as well as the options used to compile it.
    * The operating system you're running on.
    * Whether the bug occurs when running with OpenMP, MPI or both.

### Fixing a bug
* Open a pull-request (PR) with the patch and we'll look at it as soon as we get the chance (if you 
  don't know how to open a pull-request, check out 
  [this guide](https://help.github.com/articles/about-pull-requests/) from GitHub). 
* Make sure the PR description clearly describes what problem you're trying to solve and how you solved
  it. If you're fixing an existing issue, include the issue number in the description (e.g. resolves
  issue #123).
* Make sure your code compiles, passes all unit tests, and fits in with the rest of AMBiT. Importantly, 
  be consistent with the coding style in the rest of the codebase, write good comments, and follow 
  [idiomatic C++11 conventions](https://herbsutter.com/elements-of-modern-c-style/). See the
  **Development process** section of this guide for more details.

### Requesting or contributing a feature
* If you have a request for a feature, contact the developers directly with the details of your request
  and rationale behind it. Don't open a new issue, as GitHub issues are for bugs and unexpected
  behaviour. Additionally, keep in mind that we are a small team with limited resources, so we may not 
  be able to work on feature requests for a while (if at all).
* If you have a patch which introduces a new feature, email the developers with the details of your
  patch, the rationale behind it and a (small) proof of concept and we'll get back to you as soon as we
  can. Pull-requests are also welcome, but we'd prefer if you also sent us a more detailed email to go 
  with it. We work from the `dev` branch when developing new features, so please make sure your patch is
  based off of `dev`, not `master`.

Lastly, if you have any questions about contributing: ask us! We'll be happy to help get you up and 
running with developing AMBiT.

## Development process
The key thing to remember is that AMBiT follows the *Don't Be A Jerk* coding conventions: be 
considerate of other people's code, don't break anything and write useful comments. Don't Be A Jerk and 
you'll be fine.

Additionally, we try to adhere to a set of standard software engineering practices while developing 
AMBiT. While the may seem overly burdensome if you're coming from another scientific codebase, these 
practices allow us to develop solid, stable and correct software, while spending less time chasing bugs 
or untangling masses of code to integrate a new feature. Make sure you adhere to these development 
practices when writing code for AMBiT.

### Git workflow
* If you're developing a new feature, make a new Git branch for it. Branches are very cheap in Git and
  make the whole development process (from writing and testing code, to submitting pull-requests) much
  less painful. 
* Commits should be atomic (i.e. one commit for one conceptual change) but try not to make unnecessary 
  small commits, as these can clutter up the Git logs. Try to squash small, related commits together, if
  possible.
* Merging is generally less of a hassle if you regularly [rebase](https://git-scm.com/docs/git-rebase) 
  your feature branch against `dev` so it doesn't get too out-of-sync. Rebasing in 
  [interactive mode](https://git-scm.com/docs/git-rebase#_interactive_mode) has the added bonus of
  letting you squash or re-order commits to make a cleaner record of your changes.
* The [Git manual](https://git-scm.com/book/en/v2) is a good guide if you need an introduction or
  refresher on Git.

### Testing
AMBiT has a test-suite, which we require all pull-requests and commits to pass.
Breaking the build is never fun, and a good test suite lets us catch a lot of bugs well before they make
it to release.

Our testing frame work is `googletest`, which can be obtained from 
[this GitHub repository](https://github.com/google/googletest). Make sure you have this installed and
working before you start developing for AMBiT, since it's easiest to check tests as you go.

To compile the tests, do:

`scons test`

in the top-level directory. To run the tests, do:

`./ambit_test`

The test-suite will report any failing tests, which you should triage before going any further.
Additionally, make sure your feature compiles and works correctly with OpenMP, MPI, or both enabled.

### Coding conventions
* Follow the style in the rest of the code-base, including indenting with four real-spaces. Make sure to
  give your classes, variables, functions, etc useful and descriptive names. Avoid single-letter
  identifiers (like `x` or `n`) outside of loop counters, iterators or variables in lambda functions. Be
  generous (but not too verbose) with comments. Remember, code is *read* more often than it is 
  *written*, so write code accordingly.
* Don't forget to update the AMBiT user guide if you add or change a feature.
* AMBiT makes extensive use of features and idioms from C++11, as well as more general object-oriented
  design. If you're coming from a language like C or Fortran (or old C++), it's worth familiarising
  yourself with these conventions: they genuinely lead to better code. 
  [cplusplus.com](http://www.cplusplus.com/) is a good, beginner's level introduction to C++, while
  Herb Sutter's [Elements of Modern C++ style](https://herbsutter.com/elements-of-modern-c-style/) and 
  Scott Meyers' [Effective Modern C++](http://shop.oreilly.com/product/0636920033707.do) are both 
  excellent guides to modern C++ development.
* Usability is a first-class constraint for AMBiT: making the software smooth and relatively
  intuitive to use is a major key to its success. Take a moment to consider whether the feature you've
  written would be intuitive to someone who knows the physics, but has never used your code before. Try
  to eliminate "gotchas" or counter-intuitive input options where possible.

That's it! Contact the AMBiT developers if you have any questions or feedback. Otherwise: happy coding!
