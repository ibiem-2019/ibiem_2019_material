The IBIEM bootcamp and the CSP1 and CSP2 courses depend heavily on
computing so please read this information carefully. This document is
also available at
<a href="https://github.com/ibiem-2019/ibiem_2019_material/blob/master/content/misc/ibiem_computing.md" class="uri">https://github.com/ibiem-2019/ibiem_2019_material/blob/master/content/misc/ibiem_computing.md</a>

Computer
========

Many bootcamp sessions and CSP1 and CSP2 class will involve computer
work, so you *must* bring a computer with you to bootcamp, CSP1, and
CSP2. You do not need a brand-new, super-fast computer. As discussed
below, all of our computational work will be in a web-browser based
environment, so you just need to be able to run an up-to-date web
browser. If you cannot bring a computer with you, please contact us as
soon as possible so we can make arrangements for a loaner computer.

NetID and password
==================

Please be sure that you know your Duke NetID and password, you will need
it to access the Duke wireless network during bootcamp.

IBIEM computing environment
===========================

All computational work for the IBIEM will be done through a web
browser-based version of RStudio. Our RStudio instances are being run by
the Duke Office of Information Technology, so all you need on your
computer is an up-to-date version of Firefox, Chrome, or Safari.

Accessing the IBIEM computing environment
-----------------------------------------

we will be using a web based version RStudio which is hosted by Duke’s
Office of Information Technology. To access your instance:

1.  Go to
    <a href="https://vm-manage.oit.duke.edu/" class="uri">https://vm-manage.oit.duke.edu/</a>
2.  Logon with your NetID and Password
3.  Click on “Docker” in the top right.
4.  You will see a list of the available containers, find the section
    “IBIEM2019 - RStudio for bootcamp, CSP1, and CSP2 2019”
5.  Click “Click here to create your personal IBIEM2019 environment”
    -   Important: Be sure that you select IBIEM2019! There are versions
        of the IBIEM environment from previous years! You do not want
        these!!

### NCA&T Students

Sometimes the “VM-Manage” system does not work with guest NetIDs. If you
are unable to get to your IBIEM computing environment using the
instructions above, please let me know what error message you are
getting and we will work to solve the problem.

IBIEM Login confirmation
========================

In order to catch and fix any problems ahead of time, you *must* log
into the IBIEM Computing Environment as soon as possible and confirm
your successful login. You will be recieving an email later today with a
link for a “Login Confirmation” form. To complete this process please do
the following:

1.  Follow the instructions above for logging into the IBIEM Computing
    Environment

2.  Once RStudio is open in your web browser, type or paste the
    following command into the *Console* pane (on the left side of the
    RStudio window), then hit the *return* key on your keyboard:

``` r
readLines("/data/confirm.md")
```

You should see the following output:

    [1] "The login confirmation number is: #####"

where \#\#\#\#\# is a 5 digit number.

1.  You will need to type or paste the 5 digit number into the “Login
    Confirmation” form and submit it.

2.  Once you do that, you should receive a message “IBIEM Login
    Confirmed! Thank you!”
