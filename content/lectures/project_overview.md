Start Reproducibility Demo
==========================

Project Organization
====================

Class, Group, and Individual
----------------------------

> -   Class: `/data`
> -   Group: `/sharedspace/OurProject`
> -   Individual: `/home/guest`

Name that Location
------------------

> -   Code?
> -   Taxonomic Reference (e.g. Silva)?
> -   Raw Data (e.g. FASTQs)?
> -   Intermediate Files (e.g. Filtered FASTQs)?
> -   Final Results (e.g. RDS of phyloseq object)
> -   Metadata

Code
----

-   Git Repo in Home directory
-   Remote repository on Github
    -   All project members are collaborators

Raw Data
--------

-   Working Copy
    -   Group Directory
    -   Read Only
-   Archive
    -   New Data: Duke Data Service
    -   Published Data: SRA

Space: Availabile
-----------------

``` bash
df -h
```

    ## Filesystem      Size  Used Avail Use% Mounted on
    ## none            197G   61G  127G  33% /
    ## tmpfs            64M     0   64M   0% /dev
    ## tmpfs            28G     0   28G   0% /sys/fs/cgroup
    ## shm              64M     0   64M   0% /dev/shm
    ## /dev/sdd1      1008G  415G  542G  44% /home/guest
    ## /dev/sde1        50G   34G   14G  72% /tmp
    ## /dev/sdc        197G   61G  127G  33% /etc/hosts
    ## tmpfs            28G     0   28G   0% /proc/acpi
    ## tmpfs            28G     0   28G   0% /proc/scsi
    ## tmpfs            28G     0   28G   0% /sys/firmware

Space: Usage
------------

``` bash
du -h --max-depth 1 /home/guest
```

    ## 28K  /home/guest/.config
    ## 87M  /home/guest/in_class_notes
    ## 4.0K /home/guest/play
    ## 24K  /home/guest/.ssh
    ## 15M  /home/guest/.rstudio
    ## 31M  /home/guest/challenges
    ## 16K  /home/guest/R
    ## 356K /home/guest/git_demo
    ## 6.2G /home/guest/scratch
    ## 604K /home/guest/planets
    ## 4.3M /home/guest/misc
    ## 353M /home/guest/ibiem_2019_material
    ## 140K /home/guest/.cache
    ## 6.7G /home/guest

Space: Usage (continued)
------------------------

``` bash
du --max-depth 1 /home/guest | sort -nr
```

    ## 6955828  /home/guest
    ## 6454800  /home/guest/scratch
    ## 360960   /home/guest/ibiem_2019_material
    ## 88776    /home/guest/in_class_notes
    ## 31152    /home/guest/challenges
    ## 14600    /home/guest/.rstudio
    ## 4332 /home/guest/misc
    ## 604  /home/guest/planets
    ## 356  /home/guest/git_demo
    ## 140  /home/guest/.cache
    ## 28   /home/guest/.config
    ## 24   /home/guest/.ssh
    ## 16   /home/guest/R
    ## 4    /home/guest/play

Project Timeline
================

Milestones
----------

-   January 31: Group Presentations on Project Background
-   February 21: Group Presentations on Project Progress
-   April 3: Poster/Figure Critique
-   April 17: Final Presentations

Collaboration
=============

Git for Teams
-------------

-   [Team
    Conflicts](https://github.com/ibiem-2019/ibiem_2019_material/blob/master/content/lessons/bootcamp/040_git_overview.md#team-conflicts)
-   [Branching](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging)
-   [Stashing](https://git-scm.com/book/en/v2/Git-Tools-Stashing-and-Cleaning)
-   [Tagging](https://git-scm.com/book/en/v2/Git-Basics-Tagging)

Show Reproducibility Demo
=========================
