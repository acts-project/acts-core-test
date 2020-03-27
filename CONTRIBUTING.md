# Contributing to Acts

Contributions to the Acts project are very welcome and feedback on the documentation is greatly appreciated.

1. [Mailing lists](#mailing-lists)
2. [Bug reports and feature requests](#bug-reports-and-feature-requests)
3. [Make a contribution](#make-a-contribution)
    1. [Checklist for pull requests](#checklist-pull-requests)
    2. [Workflow recommendations](#workflow-recommendations)
    3. [Coding style and guidelines](#coding-style-and-guidelines)
    4. [Tips for users migrating from Gitlab] (#tips-users-gitlab)
4. [Review other contributions](#review-other-contributions)
    1. [Configuring your own PRs and PRs by people with read rights](#apply-labels)
    2. [Approving a pull request](#approving-a-pull-request)
    3. [Merging a pull request](#merging-a-pull-request)
5. [Administrator's corner](#admin-corner)
    1. [Making a new Acts release](#tag-release)

## <a name="mailing-lists">Mailing lists</a>

1. [acts-users@cern.ch](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=acts-users): Users of the Acts project should subscribe to this list as it provides
    - regular updates on the software,
    - a common place for asking any kind of questions.
1. [acts-developers@cern.ch](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=acts-developers): For discussions about
    - information about developer meetings,
    - a common place for technical discussions.

## <a name="bug-reports-and-feature-requests">Bug reports and feature requests</a>

To report an issue and before starting work, please create an issue in the [Github issue tracker](https://github.com/acts-project/acts-core/issues). A comprehensive explanation will help the development team to respond in a timely manner. Verbalizing the issue before starting work allows the other contributors to chime in and avoids disagreements how to progress.
- The title should summarize the issue
- Describe the issue in as much detail as possible in the comment

Github does not allow editing labels, assignees or setting milestone to non-members of a project with at least "Triage" permission. These will have to be set by members with Triage permission after an issue/PR is created.
Guidelines regarding labels, assignees and milestone therefore only concern members of acts-project with the necessary rights and can be ignored by others.

- Assign to yourself or leave empty
- Choose labels as appropriate
    - type of issue
    - which component is affected
    - urgency
    - fix versions
- bug reports
    - mention affected version(s)
    - issue type: "Bug"
    - a detailed description of the bug including a receipe on how to reproduce it and any hints which may help diagnosing the problem
- feature requests
    - issue type: "Improvement" or "New Feature"
    - a detailed description of the feature request including possible use cases and benefits for other users

## <a name="make-a-contribution">Make a contribution</a>

Anyone is welcome to contribute to Acts. Below is a short description how to contribute. If you have any questions, feel free to ask [acts-developers@cern](mailto:acts-developers@cern.ch) for help or guidance.

Please always fork the Acts repository you want to work on and create branches only in your own fork. Once you want to share your work, create a Pull Request (PR) (for gitlab users: equivalent to merge request) to the master branch of the upstream acts-project repository. If it is not yet ready to be merged in, please create a draft pull request (by clicking on the small arrow on the green "create pull request" button) to mark it work in progress. Once you want your branch to be merged in, request a review from the [reviewers team](https://github.com/orgs/acts-project/teams/reviewers). Once a draft merge request is reviewed, it can be merged in.

To get started with git, please refer to the [short introduction](http://git-scm.com/docs/gittutorial) as well as the [full git documentation](https://git-scm.com/doc). Tutorials as well as explanations of concepts and workflows with git can also be found on [Atlassian](http://www.atlassian.com/git/).
 
### <a name="checklist-pull-requests">Checklist for pull requests</a>
- Your branch has been rebased on the target branch and can be integrated through a fast-forward merge.
- A detailed description of the pull request is provided.
- The issue the PR closes is linked.
- All CI jobs pass.
- All newly introduced functions and classes have been documented properly with doxygen.
- Unit tests are provided for new functionalities.
- For bugfixes: a test case has been added to avoid the re-appearance of this bug in the future.
- All added cmake options were added to 'cmake/PrintOptions.cmake'.

### <a name="workflow-recommendations">Workflow recommendations</a>

In the following a few recommendations are outlined which should help you to get familiar with development process in the Acts project.

1. **Each development in its own branch of your private fork!**
Write access for developers has been disabled for developers on acts-project. Therefore, always start by creating your own fork and creating branches there. You should start a new branch for every development and all work which is logically/conceptually linked should happen in one branch. Try to keep your branches short. This helps immensly to understand the git history if you need to look at it in the future and helps reviewers of your code.
If projects are complex (e.g. large code refactoring or complex new features), you may want to use _sub_-branches from the main development branch as illustrated in the picture below.

<img src="doc/figures/sub_dev.png" alt="workflow for large feature">
1. **Never, ever directly work on any "official" branch!**
Though not strictly necessary and in the end it is up to you, it is strongly recommended that you never commit directly on a branch which tracks an "official" branch. As all branches are equal in git, the definition of "official" branch is quite subjective. In the Acts project you should not work directly on branches which are **protected** in the repository. Usually, these are the _master_ and _release-X.Y.Z_ branches. The benefit of this strategy is that you will never have problems to update your fork. Any git merge in your local repository on such an "official" branch will always be a fast-forward merge.

1. **Use atomic commits!**
Similarly to the concept of branches, each commit should reflect a self-contained change. Try to avoid overly large commits (bad examples are for instance mixing logical change with code cleanup and typo fixes).

1. **Write good commit messages!**
Well-written commit messages are key to understand your changes. There are many guidelines available on how to write proper commit logs (e.g. [here](http://alistapart.com/article/the-art-of-the-commit), [here](http://chris.beams.io/posts/git-commit/), or [here](https://wiki.openstack.org/wiki/GitCommitMessages#Information_in_commit_messages)). As a short summary:
    - Structure your commit messages into short title (max 50 characters) and longer description (max width 72 characters)!
      This is best achieved by avoiding the `commit -m` option. Instead write the commit message in an editor/git tool/IDE...
    - Describe why you did the change (git diff already tells you what has changed)!
    - Mention any side effects/implications/consquences!

1. **Prefer git pull --rebase!**
If you work with a colleague on a new development, you may want to include his latest changes. This is usually done by calling `git pull` which will synchronise your local working copy with the remote repository (which may have been updated by your colleague). By default, this action creates a merge commit if you have local commits which were not yet published to the remote repository. These merge commits are considered to contribute little information to the development process of the feature and they clutter the history (read more e.g.  [here](https://developer.atlassian.com/blog/2016/04/stop-foxtrots-now/) or [here](http://victorlin.me/posts/2013/09/30/keep-a-readable-git-history)). This problem can be avoided by using `git pull --rebase` which replays your local (un-pushed) commits on the tip of the remote branch. You can make this the default behaviour by running `git config pull.rebase true`. More about merging vs rebasing can be found [here](https://www.atlassian.com/git/tutorials/merging-vs-rebasing/).

1. **Update the documentation!**
Make sure that the documentation is still valid after your changes. Perform updates where needed and ensure integrity between the code and its documentation.

### <a name="coding-style-and-guidelines">Coding style and guidelines</a>

The Acts project uses [clang-format](http://clang.llvm.org/docs/ClangFormat.html) for formatting its source code. A `.clang-format` configuration file comes with the project and should be used to automatically format the code. There are several instructions available on how to integrate clang-format with your favourite IDE (e.g. [eclipse](https://marketplace.eclipse.org/content/cppstyle), [Xcode](https://github.com/travisjeffery/ClangFormat-Xcode), [emacs](http://clang.llvm.org/docs/ClangFormat.html#emacs-integration)). The Acts CI system will automatically apply code reformatting using the provided clang-format configuration once pull requests are opened. However, developers are strongly encouraged to use this code formatter also locally to reduce conflicts due to formatting issues.

In addition, the following conventions are used in Acts code:

- Class names start with a capital letter.
- Function names start with a lower-case letter and use camel-case.
- Names of class member variables start with `m_`.
- getter methods are called like the corresponding member variable without the prefix 'get' (e.g. `covariance()` instead of `getCovariance()`)
- setter methods use the prefix 'set' (e.g. `setCovariance(...)`)
- passing arguments to functions:
    - by value for simple data types (e.g. int, float double, bool)
    - by constant reference for required input of non-trivial type
    - by (raw) pointer for optional input of non-trivial type
    - only use smart pointers if the function called must handle ownership (very few functions actually do)
- returning results from functions:
    - newly created objects should be returned<br />
      a) as unique pointer if the object is of polymorphic type or its presence is not always ensured<br />
      b) by value if the object is of non-polymorphic type and always exists
    - existing objects (e.g. member variables) should be returned by<br />
      a) const reference for custom types with costly copy constructors<br />
      b) value in all other cases
- writing unit tests:
    - place the source code of the test in the `Tests` directory, with following options:
      a) the naming of most unit test should be `{ClassName}`+`Tests.cpp`
      b) in case specific functionality of a single class is tested in different files, use `{ClassName}`+`{TestType}`+`Tests.cpp`
      c) case the written tests concerns CPU profiling, use `{ClassName}`+`Benchmark.cpp`
    - add the test as `{ClassName}`+`UnitTest` or `{ClassName}`+`{TestType}`+`UnitTest` to the test suite for unit tests
    - add the benchmark executable as `{ClassName}`+`Benchmark` to the `bin/Profiling` folder

- Doxygen documentation:
    - Put all documentation in the header files.
    - Use `///` as block comment (instead of `/* ... */`).
    - Doxygen documentation goes in front of the documented entity (class, function, (member) variable).
    - Use the \@&lt;cmd&gt; notation.
    - Functions and classes must have the \@brief description.
    - Document all (template) parameters using \@(t)param and explain the return value for non-void functions. Mention important conditions which may affect the return value.
    - Use `@remark` to specify pre-conditions.
    - Use `@note` to provide additional information.
    - Link other related entities (e.g. functions) using `@sa`.

**Make a bugfix while working on a feature**
    During the development of a new feature you discover a bug which needs to be fixed. In order to not mix bugfix and feature development, the bugfix should happen in a different branch. The recommended procedure for handling this situation is the following:
1. Get into a clean state of your working directory on your feature branch (either by commiting open changes or by stashing them).
1. Checkout the branch the bugfix should be merged into (either _master_ or _release-X.Y.Z_) and get the most recent version.
1. Create a new branch for the bugfix.
1. Fix the bug, write a test, update documentation etc.
1. Open a pull request for the bug fix.
1. Switch back to your feature branch.
1. Merge your local bugfix branch into the feature branch. Continue your feature development.
1. Eventually, the bugfix will be merged into _master_. Then, you can rebase your feature branch on master which will remove all duplicate commits related to the bugfix.

### <a name="tips-users-gitlab">Tips for users migrating from Gitlab</a>
- The most obvious difference first: What is called Merge Request (MR) in Gitlab is called Pull Request (PR) in Github.
- Once your PR is ready to be merged, request a review by the users in the [reviewers team](https://github.com/orgs/acts-project/teams/reviewers)
- The access rights model of Github is slightly different from Gitlab. As Acts started enforcing using your own fork with the switch to Github, developers no longer have write access to the upstream repository. Unfortunately, approval of a PR by users without write access does not automatically liberate the merge button for users with write access. Instead, users with write access have to manually check if a review from an authorized person is present before merging.

## <a name="review-other-contributions">Review other contributions</a>

Acts requires that every pull request receives at least one review by a member of the reviewers team before being merged but anyone is welcome to contribute by commenting on code changes. You can help reviewing proposed contributions by going to [the "pull requests" section of the Acts (core) Github repository](https://github.com/acts-project/acts-core/pulls).

As some of the guidelines recommended here require rights granted to the reviewers team, this guide specifically addresses the people in this team. The present contribution guide should serve as a good indication of what we expect from code submissions.

### <a name="apply-labels">Configuring your own PRs and PRs by people with read rights</a>

* Check if the "request review" label is set (equivalent to WIP tag in title of a MR in Gitlab)
* Check if the "triage" label is set and configure labels, assignees and milestone for those PR
** Needs at least label "bug", "improvement", "infrastructure" or "new feature"

### <a name="approving-a-pull-request">Approving a pull request</a>

* Does its title and description reflect its contents?
* Do the automated continuous integration tests pass without problems?
* Have all the comments raised by previous reviewers been addressed?

If you are confident that a pull request is ready for integration, please make it known by clicking the "Approve pull request" button of the Gitlab interface. This notifies other members of the Acts team of your decision, and marks the pull request as ready to be merged.

### <a name="merging-a-pull-request">Merging a pull request</a>

If you have been granted write access on the Acts repository, you can merge a pull request into the Acts master branch after it has been approved.

Github may warn you that a "Fast-forward merge is not possible". This warning means that the pull request has fallen behind the current Acts master branch, and should be updated through a rebase. Please notify the pull request author in order to make sure that the latest master changes do not affect the pull request, and to have it updated as appropriate.
