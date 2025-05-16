Version numbering used in GSAS-II
=======================================================

Git commits
------------

GSAS-II versions are numbered in several ways. Every time a change to one ore more 
GSAS-II source files is committed to a repo by the git source code management system git assigns the set of changes a 40 character hexadecimal "hash" code for that version. Usually the hash code is abbreviated to the first 6-8 characters, which is still likely to be unique. We will refer to this label as the git commit number. 

Integer Tag numbers
------------------------

When a significant change is made to GSAS-II, a developer can assign an integer tag number. These numbers are intended to be assigned consecutively and thus can be used to establish when a GSAS-II version has been superseded. As of January 2025, the current tag number is 5802. Tag numbers can be looked up `in GitHub <https://github.com/AdvancedPhotonSource/GSAS-II/tags>`_.

Release version numbers
------------------------

Software packaging systems require version numbers of the form X.Y.Z, where the numbers X, Y and Z are labeled as major, minor and micro release version numbers, such as 5.1.3. X and Y are integers. Z is allowed to have a decimal place (such as 5.1.3.1) or may contain letters (such as 5.4.0rc1), but at present GSAS-II will use an integers for Z as well. The software numbering system that GSAS-II uses is intended to follow -- more or less -- the "`effective versioning <https://jacobtomlinson.dev/effver/>`_" versioning concept, as implemented below:

* The major release (starting at 5 in 2015) will be changed in GSAS-II when there are very significant changes in how the software is installed or functions that affect large portions of how the program operates. Users are advised to look carefully at the `GSAS-II home page <https://gsasii.github.io>`_ to see what has been changed.

* When relatively limited new capabilities are added to GSAS-II or a change is made in how a small section of the program functions, the minor release number is incremented, such as changing from 5.1.3 to 5.2.0. (Note that when the minor release is incremented, the micro release is set to 0.) This will usually be accompanied with a  change to the integer tag number (see above). A change in the minor release indicates that while there may be some new functionality that would be useful to learn about, users can expect the GSAS-II to function pretty much as it did before, so it is not essential to consult the documentation.

* Changes to the micro version are made when a bug is fixed or a cosmetic change has been applied and where a new distribution is intended. Users will normally not need to consider the details of an update when only the micro version number is incremented and may not want to update every time there is a new micro version released. However, if a user reports a but and has been told that the problem has been addressed, they will want the latest version. Likewise, before reporting a bug, users are asked to please make sure that error occurs with the latest GSAS-II version. The developers may become cross after tracking down something that has been reported as not working, only to find that problem was already fixed. Not every change (git commit) with GSAS-II will be assigned a new micro version. If the change is a development that is not completed or does not introduce any new functionality (for example improvements to documentation), a new git commit tag will be generated but the mini version will not be incremented. 

Tools
--------

At present three scripts are provided for developers to use to create or work with integer tag and release version numbers, to save having to review the git logs for previously assigned numbers. All of these scripts update the `git_verinfo.py` version file which is kept so that non-git installs of GSAS-II can report the currect version number.

To assign a new integer tag number (and increment the minor version number), use script:

   `.../GSASII/install/tag-version.py`

   This script is used when there are significant changes (possibly cumulative) to the
   code and a new release should be tagged both with a new integer tag number and
   minor version number. It updates the `git_verinfo.py` version file so that
   GSAS-II can quickly report the current version numbers. 

To increment only the mini version number, use script:

   `.../GSASII/install/incr-mini-version.py`

   This script is used when the have been minor changes made in GSAS-II and 
   the changes warrent an update (for example, due to a bug fix). This cannot be used
   when a new integer tag number has been assigned to the 
   current GSAS-II commit.
   It also updates the `git_verinfo.py` version file so that
   GSAS-II can quickly report the current version numbers.

   `.../GSASII/install/save-versions.py`

   This script is used to save the latest version info in the `saved_version.py` file.
   This file is created when a non-git install is made for GSAS-II so that the
   hash number for the current git version (which cannot be
   placed into the `git_verinfo.py` file) is recorded. 
   The `save-versions.py` script does not change any tags or release version numbers. 


History
--------

When the GSAS-II project was originally started, the subversion source code management system was used to track releases of code. Subversion (svn) assigns a
consecutive number every time the code was updated and this version number was used to track GSAS-II software versions. Versions up to ~5700 (early 2024) were maintained with the Argonne svn server. For much of 2024, parallel versions of GSAS-II were kept in both svn and GitHub, but updates were applied with decreasing frequency in svn, so version numbers advanced more slowly. Finally, in November, 2024 GSAS-II version 5799 was the last update made in the Argonne subversion server and update incrementing became discretionary rather than automatic each time an update was applied to the svn server.

Release version numbering starts with version 5.0.0 in early 2025. 
