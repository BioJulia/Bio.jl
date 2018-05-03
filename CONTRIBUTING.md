# Contributing to BioJulia

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to BioJulia repositories,
which are hosted in the [BioJulia Organization](https://github.com/BioJulia) on
GitHub.

These are mostly guidelines, not rules.
Use your best judgment, and feel free to propose changes to this document in a
pull request.

## Table of contents

[I don't want to read this whole thing, I just have a question!!!](#i-dont-want-to-read-this-whole-thing-i-just-have-a-question)

[What should I know about BioJulia before I get started?](#what-should-i-know-about-biojulia-before-i-get-started)
  - [BioJulia Package Maintainers](#biojulia-package-maintainers)
  - [BioJulia Administrators](#biojulia-administrators)
  - [Etiquette and conduct](#etiquette-and-conduct)
  - [Package Conventions](#package-conventions)

[How Can I Contribute?](#how-can-i-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting an Enhancement](#suggest-an-enhancement)
  - [Making Pull Requests](#pull-requests)
  - [Become a BioJulia package maintainer](#become-a-biojulia-package-maintainer)

[Styleguides](#styleguides)
  - [Git Commit Messages](#git-commit-messages)
  - [Additional julia style suggestions](#additional-julia-style-suggestions)
  - [Documentation Styleguide](#documentation-styleguide)

[Additional notes](#additional-notes)
  - [A suggested branching model](#a-suggested-branching-model)

## I don't want to read this whole thing I just have a question!!!

We understand you are excited to get involved already!
But please don't file an issue to ask a question.
You'll get faster results by using the resources below.

We have a Gitter message chat room where the community
chimes in with helpful advice if you have questions.
If you just have a question, or a problem that is not covered by this guide,
then come on over to the Gitter and we'll be happy to help.

* [Gitter, BioJulia message board](https://gitter.im/BioJulia/Bio.jl)

## What should I know about BioJulia **BEFORE** I get started?

### BioJulia Package Maintainers

In order to provide the best possible experience for new and existing users of
Julia from the life-sciences, a little bit of structure and organization is
necessary.

Each package is dedicated to introducing a specific data type or algorithm, or
dealing with a specific biological problem or pipeline.

Whilst there are some "meta-packages" such as Bio.jl, which bundle individual
packages together for convenience of installation and use, most of the BioJulia
software ecosystem is quite decentralized.

Therefore, it made sense that maintenance of the packages should also be
fairly decentralized, to achieve this, we created the role of a "Package
Maintainer".

The maintainer(s) for a given package are listed in the packages README.md file.

The maintainers of a package are responsible for the following aspects of the
package they maintain.

1. Deciding the branching model used and how branches are protected.
2. Reviewing pull requests, and issues for that package.
3. To tag releases of a package at suitable points in the lifetime of the package.
4. To be considerate and of assistance to new contributors, new community members and new maintainers.
5. To report potential incidents of antisocial to a BioJulia admin member.

**See [HERE](#additional-notes) for extra
guidance and suggestions on branching models and tagging releases.**

Package maintainers hold **admin** level access for any package(s) for which they
are listed as maintainer, and so new contributors to BioJulia should
rest assured they will not be 'giving up' any package they transfer to BioJulia,
they shall remain that package's administrator. Package maintainers also have
**push** (but not **admin**) access to all other code packages in the BioJulia
ecosystem.

This allows for a community spirit where maintainers who are dedicated primarily
to other packages may step in to help other maintainers to resolve a PR or issue.
As such, newer maintainers and researchers contributing a package to the BioJulia
ecosystem can rest assured help will always be at hand from our community.

However, if you are a maintainer stepping in to help the maintainer(s) dedicated
to another package, please respect them by first offering to step in and help,
before changing anything. Secondly, ask them before doing
advanced and potentially destructive git operations e.g forcing pushes to
branches (especially master), or re-writing history of branches.
Please defer to the judgement of the maintainers dedicated in the README of the
package.

### BioJulia Administrators

BioJulia has a select group of members in an Admin team.
This team has administrative access to all repositories in the BioJulia project.

The admin team is expected to:

1. Respond and resolve any disputes between any two BioJulia contributors.
2. Act as mentors to all other BioJulia maintainers.
3. Assist maintainers in the upkeep of packages when requested. Especially when
   more difficult re-bases and history manipulation are required.
4. Some administrators maintain the BioJulia infrastructure.
   This includes being responsible for the accounts and billing of any
   platforms used by BioJulia, and the maintenance of any hardware like
   servers owned and used by BioJulia.

### Etiquette and conduct

BioJulia outlines a [statement of etiquette and conduct](CODE_OF_CONDUCT.md)
that all members and contributors are expected to uphold. Please take the time
to read and understand this statement.

### Package conventions

First, be familiar with the official julia documentation on:

* [Packages](https://docs.julialang.org/en/stable/manual/packages/)
* [Package Development](https://docs.julialang.org/en/stable/manual/packages/#Package-Development-1)
* [Modules](https://docs.julialang.org/en/stable/manual/modules/)

Package names should be a simple and self explanatory as possible, avoiding
unneeded acronyms.

Packages introducing some key type or method/algorithm should be named
accordingly.

For example, the BioJulia package introducing biological sequence types and
functionality to process sequence data is called "BioSequences".
GitHub repository names of BioJulia packages should end in `.jl`, even though
the package name itself does not.
i.e. "BioSequences" is the name of the package, and the name of its GitHub
repository is "BioSequences.jl".

Considerate and simple naming greatly assists people in finding the kind of
package or functionality they are looking for.

File names of files containing julia code in packages should end in `.jl`.

All user facing types and functions (i.e. all types and functions
exported from the module of a package), should be documented.
Documentation regarding specific implementation details that aren't relevant
to users should be in the form of comments. Please *DO* comment liberally for
complex pieces of code!

We use [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl),
to generate user and developer documentation and host it on the web.
The source markdown files for such manuals is kept in the `docs/src/`
folder of each BioJulia package/repository.

The code in all BioJulia packages is unit tested. Such tests should be
organized into contexts, and into separate files based on module.

Files for tests for a module go into an appropriately named folder, within the
`test` folder in the git repo.

## How can I contribute?

### Reporting Bugs

Here we show you how to submit a bug report for a BioJulia repository.
If you follow the advice here, BioJulia maintainers and the community will
better understand your report :pencil:, be able to reproduce the behaviour
:computer: :computer:, and identify related problems :mag_right:.

#### Before creating a bug report:

Please do the following:

1. Check the GitHub issue list for the package that is giving you problems.

2. If you find an issue already open for your problem, add a comment to let
  everyone know that you are experiencing the same issue.

3. If no **currently open** issue already exists for your problem that has already been
   then you should create a new issue.

   > **Note:** If you find a **Closed** issue that seems like it is the same thing
   > that you're experiencing, open a new issue and include a link to the original
   > issue in the body of your new one.

#### How to create a (good) new bug report:

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/).
After you've determined [which repository](https://github.com/BioJulia)
your bug is related to, create an issue on that repository and provide the
following information by filling in [template](.github/ISSUE_TEMPLATE.md).
This template will help you to follow the guidance below.

When you are creating a bug report, please do the following:

1. **Explain the problem**

   - *Use a clear and descriptive title* for the issue to identify the problem.
   - *Describe the exact steps which reproduce the problem* in as many details as possible.
     - Which function / method exactly you used?
     - What arguments or parameters were used?
     - *Provide a specific example*. (Includes links to pastebin, gists and so on.)
       If you're providing snippets in the issue, use
       [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - *Describe the behaviour you observed after following the steps*
     - Point out what exactly is the problem with that behaviour.
     - *Explain which behaviour you expected to see instead and why.*
     - *OPTIONALLY: Include screenshots and animated GIFs* which show you
       following the described steps and clearly demonstrate the problem.
       You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on
       macOS and Windows, or [this tool](https://github.com/colinkeenan/silentcast)
       or [this tool](https://github.com/GNOME/byzanz) on Linux.

2. **Provide additional context for the problem (some of these may not always apply)**

   - *Did the problem start happening recently* (e.g. after updating to a new version)?
     - If the problem started recently, *can you reproduce the problem in older versions?*
     - Do you know the most recent package version in which the problem doesn't happen?

   - *Can you reliably reproduce the issue?* If not...
     - Provide details about how often the problem happens.
     - Provide details about under which conditions it normally happens.

   - Is the problem is related to *working with files*? If so....
     - Does the problem happen for all files and projects or only some?
     - Does the problem happen only when working with local or remote files?
     - Does the problem happen for files of a specific type, size, or encoding?
     - Is there anything else special about the files you are using?

3. **Include details about your configuration and environment**

- *Which version of the package are you using?*

- *What's the name and version of the OS you're using?*

- *Which julia packages do you have installed?*

- Are you using local configuration files to customize julia behaviour? If so...
  - Please provide the contents of those files, preferably in a
  [code block](https://help.github.com/articles/markdown-basics/#multiple-lines)
  or with a link to a [gist](https://gist.github.com/).

*Note: All of the above guidance is included in the [template](.github/ISSUE_TEMPLATE.md) for your convenience.*

### Suggest an Enhancement

This section explains how to submit an enhancement proposal for a BioJulia
package. This includes completely new features, as well as minor improvements to
existing functionality.
Following these suggestions will help maintainers and the community understand
your suggestion :pencil: and find related suggestions :mag_right:.

#### Before Submitting An Enhancement Proposal

* **Check if there's already [a package](https://github.com/BioJulia) which provides that enhancement.**

* **Determine which package the enhancement should be suggested in.**

* **Perform a cursory issue search** to see if the enhancement has already been suggested.
  * If it has not, open a new issue as per the guidance below.
  * If it has...
    * Add a comment to the existing issue instead of opening a new one.
    * If it was closed, take the time to understand why this was so (it's ok to
      ask! :) ), and consider whether anything has changed that makes the reason
      outdated. If you can think of a convincing reason to reconsider the
      enhancement, feel free to open a new issue as per the guidance below.

#### How to submit a (good) new enhancement proposal

Enhancement proposals are tracked as
[GitHub issues](https://guides.github.com/features/issues/).
After you've determined which package your enhancement proposals is related to,
create an issue on that repository and provide the following information by
filling in [template](.github/ISSUE_TEMPLATE.md).
This template will help you to follow the guidance below.

1. **Explain the enhancement**
   - *Use a clear and descriptive title* for the issue to identify the suggestion.
   - *Provide a step-by-step description of the suggested enhancement* in as many details as possible.
   - *Provide specific examples to demonstrate the steps*.
     Include copy/pasteable snippets which you use in those examples, as
     [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - If you want to change current behaviour...
     - Describe the *current* behaviour.
     - *Explain which behaviour you expected* to see instead and *why*.
     - *Will the proposed change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.

   - *OPTIONALLY: Include screenshots and animated GIFs*.
     You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on
     macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast)
     or [this tool](https://github.com/GNOME/byzanz) on Linux.

2. **Provide additional context for the enhancement**

   - *Explain why this enhancement would be useful* to most BioJulia users and
     isn't something that can or should be implemented as a separate package.

   - *Do you know of other projects where this enhancement exists?*

3. **Include details about your configuration and environment**

   - Specify which *version of the package* you're using.

   - Specify the *name and version of the OS* you're using.

*Note: All of the above guidance is included in the [template](.github/ISSUE_TEMPLATE.md) for your convenience.*

### Making Pull Requests

BioJulia packages (and all julia packages) can be developed locally.
For information on how to do this, see this section of the julia
[documentation](https://docs.julialang.org/en/stable/manual/packages/#Package-Development-1).

Before you start working on code, it is often a good idea to open an enhancement
[suggestion](#suggest-an-enhancement)

Once you decide to start working on code, the first thing you should do is make
yourself an account on [Github](https://github.com).
The chances are you already have one if you've done coding before and wanted to
make any scripts or software from a science project public.

The first step to contributing is to find the
[BioJulia repository](https://github.com/BioJulia) for the package.
Hit the 'Fork' button on the repositories page to create a forked copy of the
package for your own Github account. This is your blank slate to work on, and
will ensure your work and experiments won't hinder other users of the released
and stable package.

From there you can clone your fork of the package and work on it on your
machine using git.
Here's an example of cloning, assuming you already forked the BioJulia package "BioSequences.jl":

```sh
git clone https://github.com/<YOUR_GITHUB_USERNAME_HERE>/BioSequences.jl.git
```

Git will download or "clone" your fork and put it in a folder called
BioSequences.jl it creates in your current directory.

It is beyond the scope of this document to describe good git and github use in
more specific detail, as the folks at Git and GitHub have already done that wonderfully
on their own sites. If you have additional questions, simply ping a BioJulia
member or the [BioJulia Gitter](https://gitter.im/BioJulia/Bio.jl).

#### How to make (good) code contributions and new Pull-Requests

1. **In your code changes**

   - **Branch properly!**
     - If you are making a bug-fix, then you need to checkout your bug-fix branch
       from the last release tag.
     - If you are making a feature addition or other enhancement, checkout your
       branch from master.
     - See [here](#a-suggested-branching-model) for more information (or ask a package maintainer :smile:).

   - Follow the [julia style guide](https://docs.julialang.org/en/stable/manual/style-guide/).

   - Follow the [additional style suggestions](#additional-julia-code-style-suggestions).

   - Follow the [julia performance tips](https://docs.julialang.org/en/stable/manual/performance-tips/).

   - Update and add docstrings for new code, consistent with the [documentation styleguide](https://docs.julialang.org/en/stable/manual/documentation/).

   - Update information in the documentation located in the `docs/src/`
     folder of the package/repository if necessary.

   - Ensure that unit tests have been added which cover your code changes.

   - Ensure that you have added an entry to the `[UNRELEASED]` section of the
     manually curated `CHANGELOG.md` file for the package. Use previous entries as
     an example. Ensure the `CHANGELOG.md` is consistent with the
    recommended [changelog style](EXAMPLE_CHANGELOG.md).

   - All changes should be compatible with the latest stable version of
     Julia.

   - Please comment liberally for complex pieces of internal code to facilitate comprehension.

2. **In your pull request**

   - **Use the [pull request template](.github/PULL_REQUEST_TEMPLATE.md)**

   - *Describe* the changes in the pull request

   - Provide a *clear, simple, descriptive title*.

   - Do not include issue numbers in the PR title.

   - If you have implemented *new features* or behaviour
     - *Provide a description of the addition* in as many details as possible.
     - *Provide justification of the addition*.
     - *Provide a runnable example of use of your addition*. This lets reviewers
       and others try out the feature before it is merged or makes it's way to release.

   - If you have *changed current behaviour*...
     - *Describe the behaviour prior to you changes*
     - *Describe the behaviour after your changes* and justify why you have made the changes.
     - *Does your change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.
     - If you are implementing changes that are intended to increase performance, you
       should provide the results of a simple performance benchmark exercise
       demonstrating the improvement. Especially if the changes make code less legible.

*Note: All of the above guidance is included in the [template](.github/PULL_REQUEST_TEMPLATE.md) for your convenience.*

#### Reviews and merging

You can open a pull request early on and push changes to it until it is ready,
or you can do all your editing locally and make a pull request only when it is
finished - it is up to you.

When your pull request is ready on Github, mention one of the maintainers of the repo
in a comment e.g. `@Ward9250` and ask them to review it. You can also use Github's
review feature. They will review the code and documentation in the pull request,
and will assess it.

Your pull request will be accepted and merged if:

1. The dedicated package maintainers approve the pull request for merging.
2. The automated build system confirms that all unit tests pass without any issues.

There may be package-specific requirements or guidelines for contributors with
some of BioJulia's packages. Most of the time there will not be, the maintainers
will let you know.

It may also be that the reviewers or package maintainers will want to you to make
changes to your pull request before they will merge it. Take the time to
understand why any such request has been made, and freely discuss it with the
reviewers. Feedback you receive should be constructive and considerate
(also see [here](#etiquette-and-conduct)).

### Submitting a package to BioJulia

If you have written a package, and would like to have it listed under -
and endorsed by - the BioJulia organization, you're agreeing to the following:

1. Allowing BioJulia to have joint ownership of the package.
   This is so that the members can help you review and merge pull requests and
   other contributions, and also help you to develop new features.
   This policy ensures that you (as the package author and current maintainer)
   will have good support in maintaining your package to the highest possible
   quality.

2. Go through a joint review/decision on a suitable package name.
   This usually the original package name. However, package authors may be asked
   to rename their package to something more official and discoverable (by
   search engines and such) if it is contentious or non-standard.

To submit your package, follow these steps:

1. Introduce yourself and your package on the BioJulia Gitter channel.
2. At this point maintainers will reach out to mentor and vouch for you and your package. They will:
  1. Discuss with you a suitable name.
  2. Help you ensure the the package is up to standard, and meets the code and contribution guidelines described on this site.
  3. Add you to the BioJulia organisation if you wish to become a BioJulia maintainer.
  4. Transfer ownership of the package.

### Become a BioJulia package maintainer

You may ask the current admin or maintainers of a BioJulia package to invite you.

They will generally be willing to do so if you have done one or
more of the following to [contribute](#how-can-i-contribute) to BioJulia in the past:

1. You have [submitted a new package](#submitting-a-package-to-biojulia) to BioJulia.
2. [Reported a bug](#reporting-bugs).
3. [Suggested enhancements](#suggesting-enhancements).
4. [Made one or more pull requests](#pull-requests) implementing one or more...
    - Fixed bugs.
    - Improved performance.
    - Added new functionality.
    - Increased test coverage.
    - Improved documentation.

None of these requirements are set in stone, but we prefer you to have done one
or more of the above, as it gives good confidence that you are familiar with the
tasks and responsibilities of maintaining a package used by others, and are
willing to do so.
Any other avenue for demonstrating commitment to the community and the
GitHub organisation will also be considered.

### BioJulia members can sometimes become administrators

Members of the admin team have often been contributing to BioJulia for a long
time, and may even be founders present at the inception of the project.
In order to become an admin, one does not necessarily have to contribute large
amounts of code to the project.
Rather the decision to on-board a member to an admin position requires a history
of using and contributing to BioJulia, and a positive
interaction and involvement with the community. Any BioJulia member fulfilling
this, may offer to take on this [responsibility](#biojulia-administrators).

## Styleguides

### Git Commit messages

* Use the present tense ("Add feature" not "Added feature").
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...").
* Limit the first line to 72 characters or less.
* Reference issues and pull requests liberally after the first line.
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :green_heart: `:green_heart:` when fixing the CI build
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies
    * :exclamation: `:exclamation:` when removing warnings or depreciations

### Additional julia style suggestions

- Source code files should have the following style of header:

  ```julia
  # Title
  # =====
  #
  # Short description.
  #
  # [Long description (optional)]
  #
  # This file is a part of BioJulia. License is MIT: <link to the license file>
  ```

- Indent with 4 spaces.

- For functions that are not a single expression, it is preferred to use an explicit `return`.
  Be aware that functions in julia implicitly return the the result of the last
  expression in the function, so plain `return` should be used to indicate that
  the function returns `nothing`.

- Type names are camel case, with the first letter capitalized. E.g.
  `SomeVeryUsefulType`.

- Module names should be camel case.

- Separate logical blocks of code with one blank line. Although it is common
  and acceptable for short single-line functions to be defined together on
  consecutive lines with no blank lines between them.

- Function names, apart from constructors, are all lowercase.
  Include underscores between words only if the name would be hard
  to read without.
  E.g.  `start`, `stop`, `find_letter` `find_last_digit`.
  It is good to separate concepts in a name with a `_`.

- Generally try to keep lines below 100-columns, unless splitting a long line
  onto multiple lines makes it harder to read.

- Files that declare modules should only declare the module, and import any
  modules that it requires. Any subsequent significant code should be included
  from separate files. E.g.

```julia
module AwesomeFeatures

using IntervalsTrees, JSON

include("feature1.jl")
include("feature2.jl")

end
```

- Files that declare modules should have the same name name of the module.
  E.g the module `SomeModule` is declared under the file `SomeModule.jl`.

- When extending method definitions, define the methods with a module name prefix. E.g.

```julia
function Base.start(iter::YourType)
  ...
end

Base.done(iter::YourType, state) = ...
```

- Functions that get or set variables in a struct should not be
  prefixed with 'get' or 'set'.
  The getter should be named for the variable it gets, and the setter
  should have the same name as the getter, with the suffix `!`.
  For example, for the variable `names`:

```julia
name(node) # get node name
name!(node, "somename") # set node name
```

- When using conditional branching, if code is statement-like, an
  if-else block should be used. However if the code is expression-like
  then julia's ternary operator should be used.
  ```julia
  matches == sketchlen ? 1.0 : matches / (2 * sketchlen - matches)
  ```
  Some simple checks and expressions are also expressed using the `&&` or `||`
  operators instead of if-else syntax. For example:
  ```julia
  isvalid(foo) || throw(ArgumentError("$foo is not valid"))
  ```

## Additional Notes

### A suggested branching model

If you are a [dedicated maintainer](#biojulia-package-maintainers) on a BioJulia
package, you may be wondering which branching model to choose for development
and maintenance of your code.

If you are a contributor, knowing the branching model of a package may help
you work more smoothly with the maintainer of the package.

There are several options available, including git-flow.

Below is a recommended branching model for your repo, but it is
only a suggestion. What is best for you as the
[dedicated maintainer(s)](#biojulia-package-maintainers), is best for _you_.

The model below is a brief summary of the ['OneFlow model'](http://endoflineblog.com/oneflow-a-git-branching-model-and-workflow).
We describe it in summary here for convenience, but we recommend you check out
the blog article as a lot more justification and reasoning is presented on _why_
this model is the way it is.

#### During development

1. There is only one main branch - you can call it anything, but usually it's
   called `master`.

2. Use temporary branches for features, releases, and bug-fixes. These temporary
   branches are used as a convenience to share code with other developers and as a
   backup measure. They are always removed once the changes present on them are
   added to master.

3. Features are integrated onto the master branch primarily in a way which keeps
   the history linear and simple. A good compromise to the rebase vs. merge commit
   debate for this step is to first do an interactive rebase of the feature branch
   on master, and then do a non-fast-forward merge.
   Github now does squashed commits when merging a PR and this is fine too.

_Feature Example:_

```sh
git checkout -b feature/my-feature master

... Make commits to feature/my-feature to finish the feature ...

git rebase -i master
git checkout master
git merge --no-ff feature/my-feature
git push origin master
git branch -d feature/my-feature
```

#### :sparkles: Making new releases

1. You create a new branch for a new release. It branches off from `master` at the
   point that you decided `master` has all the necessary features. This is not
   necessarily the tip of the `master` branch.

2. From then on new work, aimed for the _next_ release, is pushed to `master` as
   always, and any necessary changes for the _current_ release are pushed to the
   release branch. Once the release is ready, you tag the top of the release branch.

3. Once the release is ready, tag the top of the release branch with a version
   number. Then do a typical merge of the release branch into `master`.
   Any changes that were made during the release will now be part of `master`.
   Delete the release branch.

_Release Example:_

```sh
git checkout -b release/2.3.0 9efc5d

... Make commits to release/2.3.0 to finish the release ...

git tag 2.3.0
git checkout master
git merge release/2.3.0
git push --tags origin master
git branch -d release/2.3.0
git push origin :release/2.3.0
```

7. Do your pushes, and go to GitHub to make your release available.

#### :bug: Hot-fixes and hot-fix releases

1. When a hot-fix is needed, create a hot-fix branch, that branches from the
   release tag that you want to apply the fix to.

2. Push the needed fixes to the hot-fix branch.

3. When the fix is ready, tag the top of the fix branch with a new release,
    merge it into master, finally delete the hot-fix branch.

_Hot-fix example:_

```sh
git checkout -b hotfix/2.3.1 2.3.0

... Add commits which fix the problem ...

git tag 2.3.1
git checkout master
git merge hotfix/2.3.1
git push --tags origin master
git branch -d hotfix/2.3.1
```

**IMPORTANT:**
There is one special case when finishing a hot-fix branch.
If a release branch has already been cut in preparation for the next release
before the hot-fix was finished, you need to merge the hot-fix branch not to
master, but to that release branch.
