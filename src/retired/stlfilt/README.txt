===========================================================================
         BD Software STL Error Message Decryptor for gcc (all versions)

 Note: Versions for MSVC 6, Visual Studio.NET 200x, Comeau C++, Intel C++,
   CodeWarrior, EDG and Borland and Digital Mars C++ are available at:
                http://www.bdsoft.com/tools/stlfilt.html

                     Copyright 2001-2008 Leor Zolman
===========================================================================

Unix/Windows:
gSTLFilt.pl (the Perl script):       Release 3.10  (01/28/2008)

Windows only:
C++.cpp (The Proxy c++ compiler):    Release 3.44  (11/25/2004)
STLTask Taskbar Icon Controller:     Release 3.07  (09/08/2004)

Written by  Leor Zolman (leor@bdsoft.com)
            BD Software
            www.bdsoft.com
                             
   BD Software offers on-site, hands-on programming seminars in
   introductory to advanced C/C++, and introductory Java, Unix 
   and Perl. You can find  more course information at the end of 
   this document, and full course descriptions at wwww.bdsoft.com.

Windows users: See QUICKSTART.txt for a brief overview and quick setup
instructions.

Unix/Linux users: The file "gfilt" is a sample shell script driver for
invoking the Decryptor. That and the script itself (gSTLFilt.pl) should
provide enough information to get you on your way.

 ****************************************************************************
 * Please do not redistribute any part of this package directly; rather,    *
 * to ensure folks get the latest version, please direct anyone interested  *
 * to download the latest distribution directly from its web page:          *
 *                                                                          *
 *      www.bdsoft.com/tools/stlfilt.html                                   *
 *                                                                          *
 ****************************************************************************

This package is updated frequently; check the web page for the date of the
latest posting, and click on the version numbers in the component table
to see reverse-chronological revision logs.

 ****************************************************************************
 *  To participate in a forum with other filter users, come visit the       *
 *  STLFilt section of the BD Software Message Board at:                    *
 *                http://bdsoft.proboards34.com/                            *
 ****************************************************************************

============================================================
What's in this archive (both Windows and Unix distributions)
============================================================

README.txt:     You're reading it: purpose, manifest, acknowledgments and 
                shameless self-promotion are all to be found here.

gfilt.BAT:      Simple command-line compiler/decryptor driver, if you'd prefer
                not to install the Proxy CL and just drive the decryptor
                manually on demand.

                Note that this batch file will NOT work under Windows95/98/ME 
                due to lack of standard *error* redirection capability in the 
                Windows9x command processor. 

                If you are running a replacement command processor for Win9x
                (such as 4DOS) that does have standard error redirection, make
                appropriate changes to the compile line and you may be able
                to get it to work.

                This error redirection issue does not arise if you use the
                "Proxy C++" command instead of the batch file.

gfilt:          Unix shell script driver for compiler/Decryptor. This script
                is smart enough to dispatch Decryptor options to the Perl
                script, and other options and arguments to the compiler.

gSTLFilt.pl:    The Perl script.

gcc-log.txt:    Revision log for C++.cpp and gSTLFilt.pl

Samples.zip/.tar:
                A slew of CPP files that draw various STL-related errors.
                Try 'em with and without the filter!

metatest.cpp:   Short sample source file illustrated the new template 
                metaprogramming-specific line wrapping features (uses Boost)

=================================
In the Windows distribution only:
=================================

QUICKSTART.txt: Fast track setup instructions, for the impatient. Contains
                just the essential instructions necessary to put the Decryptor
                into service quickly.

GLOSSARY.txt:   Definition of terms used in this software and documentation.

C++.EXE:        The C++-spoofing proxy compiler. Replaces the stock C++.EXE
                during C++ compilations, and invokes it indirectly.

STLFilt.BAT:    Low-bandwidth batch file utility to toggle the filtering
                control file. Say
                    stlfilt on
                or
                    stlfilt off
                to control filtering. Most useful when working from within
                an IDE (you can run STLFilt.BAT from a DOS window).
                
                When compiling directly from the  command line, just choose 
                between running c++2 (c++2.EXE is the recommended name for
                the copy of the native c++.exe).

                The Proxy c++'s "/NF" option may also be used to force
                "no filtering" during an individual compile.

                Note: STLFilt.BAT and STLTask.EXE cooperate, and may be run
                interchangeably to enable/disable Proxy c++ filtering.

README-STLTask.txt:
                Supplementary documentation on using the STLTask program (now
                supporting concurrent multiple compiler platforms!)

STLTask.EXE:    Taskbar-based Proxy C++ control utility.
                STLTask controls Proxy C++ installation/uninstallation, filter
                toggling and clipboard filtering for up to ten compiler
                platforms from the same instance.

                STLTask provides an alternative to STLFilt.BAT to toggle the
                filter control file. The installation/uninstallation feature
                of STLTask is provided for situations where the Proxy C++ may
                not play well with other project controlling software.
                The "Uninstall" feature makes it easy to suspend use of the
                Proxy C++ during these types of builds, if necessary, then
                bring it back for later work by selecting "Install".
                The STLTask tray icon automatically reflects the current status
                of the filtering toggle file--even when the toggle file is
                renamed using the STLFilt.BAT batch file or any other means.

Proxy-gcc.INI:  ASCII configuration file for use with C++.EXE and STLTask.EXE.
                Edit this file as appropriate for your system (remember to
                remove the leading semicolons of lines you wish to be
                active!!), place it in your Windows %SystemRoot% directory 
                C:\Windows or C:\Winnt, etc.) and its settings will override
                the hard-wired defaults built into C++.EXE and STLTask.EXE.

STLTask.INI:    Secondary (optional) configuration file for activating
                STLTask's new simultaneous multi-platform support capability.

STLFilt.BAT:    Low-bandwidth batch file utility to toggle the filtering
                control file. Say
                    stlfilt on
                or
                    stlfilt off
                to toggle whether or not the Proxy C++ performs STL message
                Decryption.
                Most useful when c++.exe is invoked automatically from within
                some kind of IDE. When compiling directly from the 
                command line, you can also simply choose between executing
                c++ (for filtering via the proxy compiler) or c++2 (to use
                the non-filtering native gcc compiler).

                The Proxy C++'s "/NF" option may also be used to force
                "no filtering" during an individual compile.

                Note: STLFilt.BAT and STLTask.EXE cooperate, and may be run
                interchangeably to enable/disable Proxy CL filtering.

C++.cpp:        Source to the Proxy C++ program.

STLTask-src.zip:
                Source code to the STLTask utility. Installation of the
                wxWindows library is required to rebuild STLTask.EXE.

STLTask-log.txt:
                Revision log for the STLTask utility.

CUSTOMIZE.txt:  Instructions for rebuilding C++.EXE and STLTask.EXE from
                their source code, and for using the STLFilt.BAT batch file
                for enabling/disabling filtering.

MinGW-install.txt
                Some notes in configuring MinGW32 gcc for use with extended
                containers.

=======
Purpose
=======

The idea is to shorten the length of STL-related error messages, and
selectively suppress several kinds of messages, so that the most vital
information out of the usual flood of diagnostics is visible front-and-center.

Lots of things you don't usually need to see are deleted, e.g.:

    --> The qualifiers "std::", "class", "struct", "__gnu_cxx", etc.
        disappear
    --> Iterators are radically shortened to either just "IT", "iter", or,
        for a container "cont", "cont::iter" (you pick which form to use).
        You can then typically deduce the type details from the remaining
        context of error messages
    --> Any functors of the form less<whatever> are deleted; others are left
        intact. Thus the default "less<...>" functors don't clutter the
        messages when dealing with associative containers
    --> strings, istreams and ostreams of <char>, their traits, etc. reduce
        to just "string", "istream", "ostream", etc. Ahhh!
    --> iostream iterators are recognized and abbreviated
    --> Allocators in type names totally disappear, and allocators in function
        parameter lists reduce to just "alloc" (good riddance...)
    --> Excessive template "candidate" lines are selectively suppressed
    --> Excessive errors from deep within the Standard Library headers are
        selectively suppressed
 
Another feature is intelligent line-wrapping for very-long type names typical
of template metaprogramming applications (if you have Boost installed, the
supplied sample source file metatest.cpp illustrates these features.)

There are certainly cases I'll have overlooked, and some I just don't know
how to filter. Feedback is welcome, but here's the disclaimer:

        The original reason I wrote this was to make it possible to teach
        C++ using the STL-first approach without scaring away students
        due to outlandishly complicated STL error messages. Thus, the filter
        addresses messages resulting from the most common coding errors (wrong
        number of arguments, wrong argument type, etc.) in standard container
        operations involving most everyday data types. If the Decryptor's
        command-line options don't provide the needed flexibility, STL power
        users are always free to disable filtering to track down errors where
        the deleted detail needs to be seen, and then re-enable filtering when
        such detail is no longer needed. Two methods are provided to enable/
        disable filtering:

            1) The STLFilt.BAT file (low-bandwidth, and calls to it can easily
               be wired into the IDE's Tools menu)
            2) the STLTask tray icon program is convenient to use, allows
                instantaneous installation/uninstallation/reinstallation of
                the Proxy C++ at any time, and provides a facility to filter
                system clipboard contents on-the-fly at any time
 

==========================================================================
           Note about the name "c++" versus the name "g++": 
In this package, I refer to the gcc compiler as "c++" since that is one of
its names. Another is "g++". I had to choose one or the other when naming
files and writing up instructions, but everything here referring to "c++"
will work just as well if you substitute "g++"... just make sure you name
files and configure parameters consistently. For example, you can rename
c++.cpp to g++.cpp (or just rename the executable to g++.exe) and
configure Proxy-gcc.INI / STLTask.INI to specify "g++" and "g++2" instead 
of "c++" and "c++2" in the appropriate sections.
==========================================================================

============
How it Works
============

1. Nuts and Bolts [Windows Only]

When installed/active, the Proxy C++ (c++.exe as distributed) is invoked
as if it were the native gcc compiler.  The Proxy C++ checks for the
existence of the controlling toggle file (FILTERING.ON). If the toggle
file is *not* detected, the Proxy C++ simply invokes the native c++
(usually renamed to C++2.EXE within gcc's bin directory) with the same
command arguments it was itself invoked with. This yields ordinary,
unfiltered error message output.

If the toggle file *is* detected AND the file type being compiled qualifies
for filtering, then the Proxy C++ sets up an interprocess pipe between the
native c++ and an invocation of Perl. The native C++'s output stream is
then piped into the standard input of the Perl process executing the
gSTLFilt.pl script, to simplify STL-related messages. The output of the Perl
script is then sent to the standard output while the process status code of
the *native* c++ process controls the subsequent behavior of any controlling 
build mechanism.


2. Rationale [Windows Only]

The Proxy C++ allows Decryptor command-line options (conveyed to the
Perl interpreter) and gcc command-line options (conveyed to the native
gcc c++ compiler) to be freely intermixed on the Proxy C++ command line. 

STLTask (the taskbar utility) has the power to install/uninstall the
filter (see the description of Active/Inactive in GLOSSARY.txt) to make it
easier to run native builds unencumbered by possible Proxy C++ bugs. This
installation/uninstallation copies either the Proxy C++ (CL.STL) or the
native c++ (C++2.EXE) over c++.EXE.  Because of the possibility of crashes
etc. in the middle of such operations, STLTask tries *real hard* to make
sure the files remain in a consistent state.  If there is any doubt as to
the current status of the filter (installed or not installed?) upon
startup, STLTask resets everything to the native gcc configuration (Proxy
c++ uninstalled) by default.


3. Configurability [Windows and Unix]

One programmer's essential details are another programmer's noise; I've put a
lot of effort into allowing users to "roll their own" feature set for the
decryption process. Be sure to carefully examine the entire "User-Configurable
Options" section near the top of the Perl script, gSTLFilt.pl, to see what
your options are. Many of those options may be controlled via the configuration
file (Proxy-gcc.INI) and command-line options to the Proxy C++, the Perl
script itself, or a driver script (see below). Some configuration options
can only be changed by editing the Perl script.

The gfilt (for Unix) and gfilt.bat (for Windows) scripts show examples of
using the options, and let you set a default set of Perl script options
via the DOPTS variable.

You may augment DOPTS variable set within the gfilt/gfilt.bat scripts by
manually setting an environment variable named GFILTOPTS to your desired
additional options (don't forget to export it under Unix).  When gfilt/
gfilt.bat run, they'll append GFILTOPTS to DOPTS.

This was done to allow quick-and- dirty option switching without having to
edit the scripts or enter a cumbersome command-line pipe sequence manually
under Windows if all you're using is the batch file.  Note that both the
Proxy C++.exe under Windows and the gfilt shell script under Unix allow
for arbitrary intermixing of compiler and Decryptor options on the command
line.


===============
The CUJ Article
===============

For background information on the origins and philosophy of the Decryptor,
see my article in the July, 2001 issue of The C/C++ Users Journal.
The article was the Web Feature for July, 2001--that means it is now
and forever available for viewing in its entirety on the CUJ web site, at
the following URL:

    http://web.archive.org/web/20041012175848/www.cuj.com/documents/s=8024/cuj0107zolman/

Note that the article describes the MSVC-specific version of the Decryptor,
so ignore installation details from the article (which are out-of-date even
for the MSVC-specific version).

        

=================================
Debugging Features [Windows Only]
=================================

C++.CPP and STLTask.EXE have debugging features built-in. C++.CPP supports
runtime debug configuration: set the DEBUG configuration option to true 
to log status messages to a file named CL-dblog.txt in the directory 
configured by the DEBUG_DIR option.

The Proxy CL option LOG_NATIVE_MSGS may be set (via the configuration file
Proxy-gcc.INI only) to force the Perl script to log its complete, UN-filtered 
input to a file named NativeLog.txt in the current directory. Sending
me a copy of this file can help me debug the Perl script when issues arise
relating to the filtering algorithm.

STLTask.cpp debugging must be configured at compile time. Set the DEBUG 
symbol to 1 and define the DEBUG_DIR symbol to the name of your
desired log file directory, and STLTask will log various status info to
a file named "debug.out" in that directory.

If you run into a situation where the Proxy C++ fails, please let
me know so I can fix it. If you understand the problem and can describe the
circumstances, that's great...if not, I'd be glad to help figure out the
problem and hopefully fix C++.cpp for the future; please reduce the project to
the minimum possible configuration that results in CL failure, zip it up
and email it to me (leor@bdsoft.com). I'll do my best to determine the 
problem and correct it. This program is as reliable as it is now *only*
because folks like you informed me about the problems (and even fixed them
for me on many occasions!)

===============
Acknowledgments
===============


Thanks to:

Scott Meyers for putting on his "Effective STL" seminar -- the event that
   inspired this project.

David Abrahams for suggesting the features relating to the /hdr:LD1, /hdr:LD2
    and /BREAK:D functionality, and then generously supporting my efforts 
    to implement those features.

David Smallberg for the Proxy compiler idea.

Thomas Becker for writing the Win32 interprocess communication code in C++.CPP

David Abrahams for contributing the long-typename-wrapping algorithm and 
assisting me in the debugging/evolving of it.


The other folks who took the time to email me with bugs and suggestions, and/or
helped with debugging:

    Thomas Becker
    Wilka
    Sam Saariste
    Scott Lewandowski
    Scott Meyers
    Tom Malcolmson
    Jan Stette
    Dominic Mathieu
    Michael Cook
    Andy Philpotts
    Rob Bishop
    John Hattan
    Derek Price
    Bill Torpey
    John Penney
    Jonathan Sambrook 
    Argos
    Matthew Douglass
    Dave Conrad
    Scott McCaskill
    Paul Suggs
    Karl Lean
    Benoit Goudreault-Emond
    Michael Griggs
    Marco Welti
    Ed Kaulakis
    Shing-Fat Fred Ma
    George Katsirelos
    Graeme Prentice
    Larry Evans

=======
Courses
=======

"Programming is not a commodity" (at least not in my book).

To promote excellence in software development technology, BD Software
offers the following on-site, hands-on seminars for programmers:


==> C++ for non-C Programmers
                   Very intensive, 5-day workshop teaching C++ for those who
                   have never programmed in C. Seminar by me, based on
                   the Stephen Prata textbook "C++ Primer Plus" (4th edition).
                   Prerequisite: *SOME* programming experience, but it
                        doesn't matter in what language.

==> An Effective Introduction to the Standard Template Library (STL)
                   Course by Scott Meyers, who needs no introduction.
                   Teaches the STL to experienced C++ programmers. 4 Days.
                   Prerequisite: Real-world programming experience in C++. 

==> C++ and Object-Oriented Programming
                   5-Day course by Dan Saks, another C++ guru you probably
                   know about.
                   Teaches C++ assuming GOOD (working) knowledge of C. More in
                   depth on OO and C++ subtleties than the first course above.

==> Advanced C Programming
                  This five day, hands-on C training is designed to bring
                  practicing  C programmers up to the next level of C
                  expertise. Since one area where C syntax and semantics
                  present a major hurdle is in the understanding 
                  of declaration syntax, the course leads participants in the
                  incremental design and implementation of their own limited C 
                  declaration parser, written (of course) in C itself. After
                  going through this process, you'll most likely understand
                  declarations almost as well as Dan Saks does (which is
                  saying a lot).

==> Introduction to Programming in C
                   My own C programming course. Some programming experience
                   (not necessarily in C) is the only prerequisite.
                   4 or 5 days; 5 preferred.
                   What can I say? Not too many folks are still learning
                   "just" C anymore, but if you have two weeks available to
                   learn C++ in, it works out a lot better to take C first and
                   then the C++ and OOP course above than it does to cram
                   both C and C++ into one week. That is, at least until
                   I write a course based on the STL-first approach...

==> Object-Oriented Programming in Java
                   My (I'm co-author) Java course. Does NOT presume any C/C++
                   experience. Uses all public domain / shareware tools (Sun's
                   JDK, ModelWorks' JPadPro (awesome IDE!). 5 Days.
                   Everyone gets a CD full of Java goodies.
                   Again, presumes only some programming experience in SOME
                   language.

==> Introduction to Unix (2-3 days)
==> Korn Shell Programming (3 days)
==> Perl Programming (3 days)

                   All three Unix courses above written by Danette Morris.
                   The Perl course is mostly platform-independent, but
                   includes some Unix-specific functionality.


To give feedback on this package, or get additional information on having one
of the courses listed above delivered at your site, please contact:

    Leor Zolman
    BD Software
    74 Marblehead Street
    North Reading, MA 01864-1527
    voice/FAX: 978-664-4178
    email: leor@bdsoft.com
    web: www.bdsoft.com

===============================================================================
