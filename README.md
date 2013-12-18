idhoa
=====

Software for Decoding of High Order Ambisonics to Irregular Layouts


## What is IDHOA? ##


IDHOA it is a small program that generates the Ambisonics decoding coefficients for irregular layouts up to fifth order. 
(At the moment not for 2D/planar layouts, only for 3D ones)

It is made of four python files with different purposes:
-- constants.py - contains all the initialization parameters and constants
-- functions.py - contains all the accessory functions 
-- plotting.py - contains only the plotting functions
-- main.py - is the main program that uses constants from constants.py and calls functions from functions.py and plotting.py



## Warning ##

This is not an out-of-the-box tool. A good understanding of Ambisonics decoding is required. Moreover, some tweaking of the parameters and a good understanding of the plots is required to get a "good" result. 

Some of the things that the program does are not necesarily correct or wrong: it depends on the USE that YOU make.


## Dependancies ##

The program depends on an external minimization library.
In principle you can use your favorite minimization library, after proper adaptation of the code. 
I tried several minimization libraries and I ended up using nlopt because it is fast enough and has a well documented api for python. [NLopt python reference.](http://ab-initio.mit.edu/wiki/index.php/NLopt_Python_Reference)  (note: in my code the implementation of constraints is still missing)
It might be that in the future I will move to something like minuit (faster? better handling of constraints?).

### How to install nlopt ###
[NLopt wiki: Installation](http://ab-initio.mit.edu/wiki/index.php/NLopt_Installation)
"By default, NLopt compiles as a static library. This means that each program you link to NLopt will make a separate copy of the library, which wastes a little disk space. The alternative is to compile NLopt as a shared library (also called a dynamic-link library)."

    ./configure --enable-shared

and then 

    make
    sudo make install

as usual.

### Python dependancies ###
The code depends on three python libraries: scipy, numpy, matplotlib.

How to install pylab, numpy, matplotlib and configure matplotlib
http://www.scipy.org/install.html
In Ubuntu and Debian:

    sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose 



## How to use IDHOA ##
1. First you have to modify constants.py to fit your needs.

* Information that you have to provide
The angular coordinates (theta and phi) of the loudspeakers, in the usual acoustics reference frame.


    * Parameters
        * Basic parameters
            1. PHI - list of phi coordinates of the loudspeakers

            2. THETA - list of theta coordinates of the louspeakers (in the same order of phi coords)

            3. DEC - available decoding schemes: basic, maxRe, phase.

            * The basic decoding initialize the minimization algorithm with the "naive" of "pinv" coefficients.
            If you selected the basic decoding then you probably don't want the minimization algorithm to search for a new minimum... and maybe you would prefere to get the pinv coefficients. Just modify the main.py accordingly.
            That said you probably would prefere to use "maxRe" or "phase" optimizations.

            * maxRe - implements the "max r_e" criterium. Tries to maximize the modulus of the energy vector r_e. (Look appendix A.4.2 Jérôme Daniel PhD thesis)

            * phase - implements the "in-phase" criterium, where all the loudspeakers are in-phase (i.e. the gains are >=0). (Look appendix A.4.3 Jérôme Daniel PhD thesis)


            4. DEG - the degree of ambisonics decoding you want

            5. SEED - gives the initial number to the function that generates the evaluation point over the sphere. A reasonable number is between 14 and 20. (You get between 245 and 503 points over the sphere respectively) 


    * Flags
        1. MUTESPKRS - tries to put to zero the smallest coefficients (at the moment if smaller than 3*10e-4). Pros: This can help to mute some misplaced spekers that are not usefully contributing to the decoding. Cons: time consuming because it runs serveral times the minimization algorithm. (TBD: write a note on the exit strategy from the loop) 

        2. AUTOREM - removes the sampling points (over the sphere) that are too far from the loudspeakers. Use this flag ONLY if you layout is not fully spherical but has some areas where there are no louspeakers. Pros: The algorithm doesn't try to optimize in regions where there are no loudspeakers. This is good because otherwise the existing loudspeakers will have to try to compensate for the missing ones. Cons: ?. It's a binary mask.

        3. PREFHOM - means "Prefere Homogeneous", if true tries to smooth the differences between the coefficients. Pros: Cons:

        4. WFRONT - gives priority to the front in optimizing the objective function

        5. WPLANE - gives priority to the horizontal plane in optimizing the objective function

        6. WBIN - you can implement a "binary mask" to remove unused areas of the sphere. This in principle can be done automatically just by using AUTOREM

        7. WAUTOREM - performs almost the same operation that AUOTREM does but using weights into the objective function. This way -if you modify the weight coefficients in "function" inside functions.py- you can adjust the weight that you give to the missing area (does this have any advantage? maybe). If you just want to remove the uncovered area then it's faster to use AUTOREM.



2. Then you might want to have a look at the most important parts of the main.py
The most important part is the MINIMIZATION section.
* a "while" loop repeats the minimization process until some condition is met (this if AUTOREM flag is true). Why? Because AUTOREM tries to set to zero small coefficients after the minimization process. Then we perform the minimization another time, to be sure that the other coefficients adjust to the new added constraint ... and maybe some other coefficient will become small enough to be neglected.
* inside the loop the real minimization is performed. 
    * First the nlopt algorithm is chosen. Several choices are possible ([have a look at this](http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms)). I choose the LN\_SBPLX because it is both fast and accurate for the purpose. Note that you might want to change the exit condition of the a) while loop accordingly to the characteristics of your algorithm.
    Then you ask nlopt to minimize the objective function called "function" (defined in functions.py), and some parameters are set (just follow the [instructions](http://ab-initio.mit.edu/wiki/index.php/NLopt_Python_Reference))
    * The (ambisonics) coefficient values are bounded to [-1,1] interval. And are put to zero if they are smaller than 3\*10e-4. This can happen if they are on a node of a SH or if the result of the minimization (\>=second loop) gives a "close to zero" value.
    * "exit condition" This condition is the one that breaks the loop. It might depend on the type of algorithm that you use (not all of them are completely deterministic). For the LN\_SBPLX I choose to exit the program when the value of the objective function is exactly the same as the previous loop. Or if the value of the objective function increased more than 1 in the last loop (then you might want to carefully check what is happening to your particular minimization, something wrong is happening). (Note: indeed it is not impossible at all that the value of the objective function increases during the process of putting to zero some coefficients)



3. The objective function.
The objective function is defined in functions.py and it is called "function". [...] I'll write more on this, but for the moment you can have a look at [this](http://www.aes.org/e-lib/browse.cfm?elib=16818)




## Utility generatingambdecpreset.pl ## 

This utility generates an ambdec preset given an header file and two coefficient files.
The header contains the usual preambleof ambdec configuration files, /description, /version, /dec\*, /opt\* and /speakers/.
The coefficients files are the output of IDHOA. The first file is used for lowfrequencies and the second is used for high frequency matrix.



## Troubleshooting ##

If, after installing nlopt, the nlopt library is not found you have to execute (once) 

    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib/

before running the program main.py



## Notes ##
### Why AUTOREM makes sense? ###
For several reasons. One interesting reasons is that it may happen that your layout is not only irregular but also "strange"... it is not fully spherical and it is not exactly an half sphere (a dome). So how do you calculate how many orders of Ambisonics you can aim to reproduce? 
This can be a question that can have a difficult answer. If you just run the program on your "strange" layout it will calculate the coefficients up to an order in which the coefficients are "sufficiently" different from zero. By looking at the coefficients values you can identify up to which order they are different from zero and then decode for that order. (btw I suggest to re-run the code specifying the lower order)
