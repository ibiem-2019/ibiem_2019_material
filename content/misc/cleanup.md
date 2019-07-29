Clean up /tmp
=============

How full is `/tmp`? You can check with the following command. If you see 100% (or something close to it, like 99%) under `Use%`, then `/tmp` is full.

``` bash
df -h /tmp
```

    ## Filesystem      Size  Used Avail Use% Mounted on
    ## /dev/sdd1        50G  2.4G   45G   6% /tmp

If it is full, please go through the following steps to clean up any large files that you no longer need. Note that /tmp is in a space that is shared by everyone, so it might not be your fault. If you clean up your `/tmp` and it doesn't significantly decrease `Use%`, then let the TAs and/or instructors know so that we can initiate a class-wide cleanup.

The following command shows you all the subdirectories in your `/tmp`, sorted by size from largest to smallest. The only subdirectories here are from R (starting with `Rtmp`), Rstudio (`rstudio-rsession`), and tmux (`tmux-1000`). We don't want to mess with any of these directories, and anyway, they aren't that big - the sizes are listed in kb. The only large directory is `/tmp` itself (102432KB is about 100MB), which means that there is one or more large files directly in `/tmp`.

``` bash
du --all --max-depth 1 /tmp | sort -nr
```

    ## 183140   /tmp
    ## 102400   /tmp/big_file.txt
    ## 38900    /tmp/scratch
    ## 32232    /tmp/emp-paired-end-sequences
    ## 9540 /tmp/casava-18-paired-end-demultiplexed.zip
    ## 16   /tmp/history.txt
    ## 8    /tmp/rstudio-rsession
    ## 4    /tmp/tmux-1000
    ## 4    /tmp/RtmprHRMiE
    ## 4    /tmp/Rtmpp7DGJR
    ## 4    /tmp/RtmpnOONaM
    ## 4    /tmp/RtmpKjj6re
    ## 4    /tmp/RtmpkaZhYt
    ## 4    /tmp/RtmphyN8t6
    ## 4    /tmp/RtmpHFHn8v
    ## 4    /tmp/Rtmpe9BhCZ
    ## 4    /tmp/Rtmp6gC0ey

The command we ran above is useful for finding out which directories are using the most storage, but since sizes are listed in kb, its harder to think about with larger files. We can run `du` with the `--human-readable` option to convert file sizes to MB and GB, where appropriate (the downside is that we can't sort to easily find the largest). This tells us that `/tmp/big_file.txt` is 100MB, so we probably want to get rid of it.

``` bash
du --all --max-depth 1 --human-readable /tmp
```

    ## 4.0K /tmp/RtmpKjj6re
    ## 4.0K /tmp/RtmpkaZhYt
    ## 4.0K /tmp/RtmprHRMiE
    ## 16K  /tmp/history.txt
    ## 8.0K /tmp/rstudio-rsession
    ## 38M  /tmp/scratch
    ## 100M /tmp/big_file.txt
    ## 9.4M /tmp/casava-18-paired-end-demultiplexed.zip
    ## 4.0K /tmp/Rtmpp7DGJR
    ## 4.0K /tmp/RtmpnOONaM
    ## 4.0K /tmp/tmux-1000
    ## 4.0K /tmp/RtmpHFHn8v
    ## 4.0K /tmp/Rtmpe9BhCZ
    ## 4.0K /tmp/Rtmp6gC0ey
    ## 4.0K /tmp/RtmphyN8t6
    ## 32M  /tmp/emp-paired-end-sequences
    ## 179M /tmp

Another tool is `ls`. It has the advantage that it can give us human readable file sizes, sorted by size. The disadvantage is that it doesn't give us acurate directory size information. The following command shows you all the files in your `/tmp`, sorted by size from largest to smallest. We want to get rid of any file we don't need, but particularly that one at the top that is 100M. &gt; If you are uncertain what can safely be deleted, feel free to check with a TAs and/or instructor.

``` bash
ls -lShap /tmp
```

    ## total 110M
    ## -rw-r--r--   1 guest users 100M Feb 19 05:37 big_file.txt
    ## -rw-r--r--   1 guest users 9.4M Oct  3 15:23 casava-18-paired-end-demultiplexed.zip
    ## -rw-r--r--   1 guest users  13K Feb 11 21:47 history.txt
    ## drwxrwxrwt  15 root  root  4.0K Feb 19 05:37 ./
    ## drwxr-xr-x 176 root  root  4.0K Feb 18 15:27 ../
    ## drwxr-xr-x   2 guest users 4.0K Jan 15 19:40 emp-paired-end-sequences/
    ## drwxrwxrwt   2 root  root  4.0K Feb 19 05:00 rstudio-rsession/
    ## drwx------   2 guest users 4.0K Feb 19 05:31 Rtmp6gC0ey/
    ## drwx------   2 guest users 4.0K Nov 16 06:45 Rtmpe9BhCZ/
    ## drwx------   2 guest users 4.0K Jan 14 15:43 RtmpHFHn8v/
    ## drwx------   2 guest users 4.0K Feb 19 05:37 RtmphyN8t6/
    ## drwx------   2 guest users 4.0K Feb 15 03:27 RtmpkaZhYt/
    ## drwx------   2 guest users 4.0K Jan 25 04:39 RtmpKjj6re/
    ## drwx------   2 guest users 4.0K Jan 20 21:23 RtmpnOONaM/
    ## drwx------   2 guest users 4.0K Feb 11 17:39 Rtmpp7DGJR/
    ## drwx------   2 guest users 4.0K Nov 15 22:56 RtmprHRMiE/
    ## drwxr-xr-x   6 guest users 4.0K Feb  8 19:59 scratch/
    ## drwx------   2 guest users 4.0K Feb 11 21:38 tmux-1000/

The following will delete `/tmp/big_file.txt`. We can also do it from the Rstudio file pane by checking the box and clicking "Delete" in the "Files" pane

``` bash
rm -f /tmp/big_file.txt
```

Clean up /home/guest
====================

We also want to check our home directory. If you run the following chunk and see 100% (or something close to it, like 99%) under `Use%`, then our space is full. Please clean up any large files that you no longer need.

``` bash
df -h /home/guest
```

    ## Filesystem      Size  Used Avail Use% Mounted on
    ## /dev/sde1      1008G  779G  179G  82% /home/guest

As with `/tmp`, we can use `du` and `ls` to find large files.

``` bash
ls -lShap ~/ | grep -v / 
```

    ## total 128K
    ## -rw-r--r--  1 guest users  17K Feb 18 19:42 .Rhistory
    ## -rw-------  1 guest users 9.8K Feb 18 23:24 .bash_history
    ## -rw-------  1 guest users  756 Jan 25 06:17 .viminfo
    ## -rw-r--r--  1 guest users  182 Jan 20 21:32 .wget-hsts
    ## -rw-r--r--  1 guest users   50 Oct 12 02:59 .gitconfig
    ## -rw-------  1 guest users   35 Jan 15 16:55 .lesshst

It doesn't look like I have any particularly large files in my home directory.

``` bash
du --all --max-depth 1 ~/ | sort -nr
```

    ## 959044   /home/guest/
    ## 726880   /home/guest/IBIEM_2018_2019_https
    ## 165740   /home/guest/solutions
    ## 32044    /home/guest/IBIEM_2018_2019
    ## 22336    /home/guest/scratch
    ## 4916 /home/guest/IBIEM_2017_2018
    ## 2616 /home/guest/.rstudio
    ## 1652 /home/guest/ncbi
    ## 528  /home/guest/challenge_7
    ## 508  /home/guest/challenge_6
    ## 428  /home/guest/challenge_8
    ## 400  /home/guest/challenge_5
    ## 364  /home/guest/.cache
    ## 260  /home/guest/challenge_4
    ## 256  /home/guest/challenge_2
    ## 28   /home/guest/.config
    ## 20   /home/guest/.ssh
    ## 20   /home/guest/.Rhistory
    ## 12   /home/guest/R
    ## 12   /home/guest/.bash_history
    ## 4    /home/guest/.wget-hsts
    ## 4    /home/guest/.viminfo
    ## 4    /home/guest/tmp
    ## 4    /home/guest/.lesshst
    ## 4    /home/guest/.gitconfig

It looks like I should do some clean up of subdirectories! Generally it is a good idea to start with the biggest subdirectory, using `du` and `ls` to find out what is taking the most space and cleaning that up first, for example, we should first look at `~/IBIEM_2018_2019_https` :

``` bash
du --all --max-depth 1 ~/IBIEM_2018_2019_https | sort -nr
```

    ## 726880   /home/guest/IBIEM_2018_2019_https
    ## 703600   /home/guest/IBIEM_2018_2019_https/.Rproj.user
    ## 12124    /home/guest/IBIEM_2018_2019_https/.git
    ## 11116    /home/guest/IBIEM_2018_2019_https/content
    ## 16   /home/guest/IBIEM_2018_2019_https/.Rhistory
    ## 4    /home/guest/IBIEM_2018_2019_https/README.Rmd
    ## 4    /home/guest/IBIEM_2018_2019_https/README.md
    ## 4    /home/guest/IBIEM_2018_2019_https/IBIEM_2018_2019_https.Rproj
    ## 4    /home/guest/IBIEM_2018_2019_https/.gitlab-ci.yml
    ## 4    /home/guest/IBIEM_2018_2019_https/.gitignore
    ## 0    /home/guest/IBIEM_2018_2019_https/misc
