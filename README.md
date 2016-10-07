    #  _________. ______. __________________.__
    # <\______   \\      \\____ _ \.___ _   |  |__ ____.    _____. .__.
    #   |  _/| ._//   ^   \|     _/|     ___|  |  \\__  \  /  ____/ __ \
    #   |   \|   /    |    |>   |  \    |   |   Y  \/ __ \_\___ \\  ___/
    #   |______  \____|__  |\___/ . \___/ . |___|  (____  /____  >\___  >
    #     .    \/ .  .   \/ .   .  .     .    .  \/ .   \/     \/ .   \/  v1.3
    #               .   .     .             .               .     . .
    # Copyright (c) 2016, Lloyd T. Elliott and Yee Whye Teh.    .
    #                     .                     .                .  .
    # All rights reserved.             .


# BNPPhase:
- Impute genetic data using a nonhomogeneous and nonparametric HMM.

- Genetic sequence data are well described by hidden Markov models (HMMs) in
 which latent states correspond to clusters of similar mutation patterns.
 Theory from statistical genetics suggests that these HMMs are nonhomogeneous
 (their transition probabilities vary along the chromosome) and have large
 support for self transitions. BNPPhase is a new nonparametric model of
 genetic sequence data, based on the hierarchical Dirichlet process, which
 supports these self transitions and nonhomogeneity. BNPPhase provides a
 parameterization of the genetic process that is more parsimonious than other
 more general nonparametric models which have previously been applied to
 population genetics. BNPPhase uses truncation-free MCMC inference derived
 through a new auxiliary sampling scheme for Bayesian nonparametric HMMs.
 The popular finite model fastPHASE can be seen as a parametric truncation of
 the BNPPhase model.

- Detail about the BNPPhase model is given in &lsquo;L.T. Elliott, Y.W. Teh. A
 nonparametric HMM for genetic imputation and coalescent inference. 2016.
 Electronic Journal of Statistics. 10(3).&rsquo;

- BNPPhase is released under the BSD 2-clause license, which is provided in
 ``LICENSE.md``. The BNPPhase executable is bundled with a Scala library and
 libraries from the Apache Commons. Scala is provided under the BSD 3-clause
 license (provided in ``classes/scala/LICENSE.md`` and the libraries from the
 Apache commons are provided under the Apache Software License 2.0 (provided in
 ``classes/commons/LICENSE.md``. 

- ``perl`` and ``java`` are required in order to run BNPPhase. This software has
 been tested with ``perl`` version 5.20 and ``java`` openjdk version 1.8, on a
 Linux system running Fedora 22. This software should support a wide variety of
 systems (including cygwin) and it should support many earlier or later versions
 of perl and java. If you encounter problems running BNPPhase please contact us
 using the information at the bottom of this manual and we would be happy to
 assist you.


--------------------------------------------------------------------------------

## USAGE:
- The following command performs imputation using the BNPPhase model. Advanced
 usage (in which control is given over all aspects of the prior, and also over
 which statistics are collected from the MCMC samples) is provided in the
 ``ADVANCED USAGE`` section below.

    ./BNPPhase -input <IN-FILE> -output <OUT-DIR> -impute [-iters <INT>]
      [-burnin <INT>] [-init <INT>] [-thin <INT>] [-random_restarts <INT>]
      [-seed <INT>]

- Performs imputation on the phased genotype file IN-FILE and stores the
 result in the directory OUT-DIR (this directory is created if it does not
 already exist). Imputed values are found by averaging over collected MCMC
 iterations.

 ``-impute`` This flag specifies that BNPPhase should perform imputation on
 missing genotypes.
 
 ``-iters <INT>`` Set the number of MCMC iterations to collect to ``INT``.
 Default value is ``60``.
 
 ``-burnin <INT>`` Set the number of MCMC iterations to perform before starting
 collection of MCMC iterations to ``INT``. Default value is ``200``.
 
 ``-init <INT>`` Set the number of MCMC iterations to do after the SMC
 initialization in which the trajectories are fixed and the parameters are
 resampled to ``INT``. These MCMC iterations occur before the burnin
 iterations. Default value is ``10``.
 
 ``-thin <INT>`` Set the number of MCMC iterations to perform and discard
 between each collected iteration to INT. Default value is ``1``. Note that the
 total number of MCMC iterations performed in which the trajectories are
 resampled is ``burnin + (thin + 1) * iters``.
 
 ``-random_restarts <INT>`` The default number of random restarts is ``1``.
 
 ``-seed <INT>`` Set the seed for the random number generator to ``INT``. If
 this flag is not provided, then a seed will be read from the random device
 ``/dev/random``.


## INPUT FORMAT:

- The input file must be a matrix of binary ``0/1`` coded SNPs. Each column
 represents a SNP, and each row represents a chromosome of a subject. The
 input is assumed to be phased: so if a diploid chromosome is being considered,
 then each subject must occupy two rows of the input file, and each row must be
 a ``0/1`` sequence representing the alleles for each of the two phased
 chromosomes for that subject. The input file must be space separated.

 Unobserved or held out alleles must be represented by a question mark ``?``.
 An example file is provided as follows:

- ``example/input.txt``

          0 0 0 1 0
          0 0 0 1 0
          0 0 0 1 ?
          1 0 0 1 0
          1 0 0 1 ?
          1 0 ? 1 ?
          1 0 ? 0 1
          1 0 ? 0 1
          1 0 1 0 1
          1 0 1 0 1

 In this example, ``10`` sequences are provided and are marked at 5 locations.
 The first two sequences are fully observed, and the third sequence has the
 ``5``-th SNP held out.


## OUTPUT FORMAT:

- The command specified in the ``USAGE`` section above creates a directory
 called ``OUT-DIR`` with the following structure:

          OUT-DIR/cmd
          OUT-DIR/solution/000000

 The file ``OUT-DIR/cmd`` contains a copy of the command line arguments used to
 specify the call to the BNPPhase software. The file ``OUT-DIR/solution/000000``
 contains calls for the alleles from the input file which are either held out or
 unobserved. The columns and rows of ``OUT-DIR/solution/000000`` match the
 columns and rows of the input file. Alleles that are observed in the input file
 are replaced with a ``-1`` in the file ``OUT-DIR/solution/000000``.


## EXAMPLE:

- The following command performs imputation on the example file provided with
 this software:

          ./BNPPhase -input example/input.txt -output example/output -impute 

 After running this command, the accuracy of BNPPhase's imputed solution can be
 computed:

          ./bin/accuracy example/validate.txt example/output/solution/000000

 The ``bin/accuracy`` command compares the calls of BNPPhase to a specification
 of the true values of alleles held out from ``example/input.txt`` (in
 ``example/validate.txt``), and outputs the proportion of calls that are
 correct. (If ``bin/accuracy`` detects an error in the alignment or
 specification of the files passed to it, then it outputs a ``?`` or a ``X``).
 The format of the validation file is the same of the format specified in the
 ``INPUT FORMAT`` section. All observed alleles must be marked as unobserved
 (i.e., with an ``?``) in the validation file. The above accuracy command
 should produce output which looks something like this:

          0.83333333

 This indicates that the accuracy of BNPPhase on the example dataset is
 ``~83%``.


## ADVANCED USAGE:

    ./BNPPhase

      # Basic options
        -input <IN-FILE>
        -output <OUT-DIR>
        [-impute]
        [-iters <INT>]
        [-burnin INT>]
        [-init <INT>]
        [-thin <INT>]
        [-random_restarts <INT>]
        [-seed <INT>]

      # System and I/O options
        [-m <INT>(G|M)]
        [-a]

      # MCMC options
        [-calls]
        [-statistics]
        [-mosaic]
        [-trace]

      # Model parameterization
        [-mass (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)]
        [-varmass (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)]
        [-mu (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)]
        [-varmu (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)]
        [-gamma (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)]
        [-r (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)]
        [-b (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)]
        [-likelihood <STRING>]
        [-mass_prior <STRING>] 
        [-mu_prior <STRNIG>]
        [-fixmass]
        [-fixmu]
        [-fixgamma]
        [-fixb]
        [-fixr]

- MCMC options. The following flags control what summary statistics are saved to
 disk in collected MCMC iterations.

 ``-trace`` If this flag is specified, calls for every iteration of the MCMC
 will be recorded and saved to the ``OUT-DIR/trace/`` directory. Here,
 ``<OUT-DIR>`` is the directory provided to the ``-output`` flag. A file for
 each iteration of the MCMC (including burnin and thinned samples) are
 reported in this directory. The probabilities used to produce these calls are
 given by the MCMC states (i.e., MCMC states are not averaged over multiple
 iterations to produce the calls in the ``OUT-DIR/trace/`` directory). The
 initial MCMC iterations (in which the parameters are updated with fixed
 trajectories) are not included in the ``OUT-DIR/trace/`` directory. The format
 used for the trace files is the same as the format described in the
 ``OUTPUT FORMAT`` section (i.e., with ``-1`` for sample/location pairs that
 appear in the input file).

 ``-mosaic`` If this flag is specified, the cluster numbers of the top level
 Dirichlet process to which each sequence is assigned at each location are
 recorded during the MCMC, forming a mosaic pattern of haplotypes [2]. This
 mosaic pattern is saved to the directory ``<OUT-DIR>/mosaic/``.

 ``-statistics`` If this flag is specified, BNPPhase will record additional
 statistics for each MCMC iteration. These statistics are saved in the 
 following files in the ``<OUT-DIR>`` directory: ``beta``, ``gamma``, ``r``,
 and ``statistics``. The ``beta`` file contains the value of the location
 dependent mean minor allele frequency. This is the parameter ``b`` from [1].
 Each line contains a space separated list of the values, one for each marker.
 Each MCMC iteration is described by one of the lines (including the initial
 state, which is on the first line). The ``gamma`` file contains the location
 dependent minor allele frequency pseudo-observations, in the same format as
 the ``beta`` file. These two files are only meaningfull if the ``-likelihood``
 flag does not specify an allele likelihood other than the default
 beta/Bernoulli likelihood (the ``-likelihood`` flag is described below).

 The ``r`` file contains the values of the jump rates at each iteration (with
 one line per iteration). The number of jump rates is one less than the number
 of markers (there is no jump rate associated with the first marker), and that
 is reflected in the number of columns of ``r``.

 The sampled values of the latent parameters and some marginal statistics for
 the BNPPhase model are contained in the ``statistics`` file. Again, there is
 one line per iteration and in addition the first line is a header describing
 the columns.

 ``-calls`` If this flag is specified, then for each MCMC iteration, calls are
 produced in the ``<OUT-DIR>/calls/`` directory by combining marginals for
 each previously collected MCMC iteration. The format of the files in this
 directory is the same as the format described in the ``OUTPUT FORMAT``
 section.

- System and I/O options:

 ``-m <INT>(G|M)`` This flag sets the internal memory used by BNPPhase to
 ``INT``. The suffix ``G`` or ``M`` denotes ``GBs`` or ``MBs`` respectively. By
 default, the amount of memory used by BNPPhase is ``1G``.
 
 ``-a`` This flag spcifies that BNPPhase should conduct internal assertions and
 verifications. If BNPPhase crashes or produces unexpected results, please
 rerun with the ``-a`` flag and send us a bug report along with a stack trace
 (if any).

- Model specification. The following flags control the initialization and range
 of latent variables used by the BNPPhase model. There are three ways of
 specifying each flag: ``<INIT>`` (without angle brackets) will set the initial
 value of the variable to ``INIT`` (a real number). ``<MIN>:<MAX>`` will set the
 minimum and maximum values for a variable to ``MIN`` and ``MAX`` respectively.
 The prior on the variable will be truncated to the range ``[MIN, MAX]``.
 ``<MIN>:<INIT>:<MAX>`` will set both the initial value of the variable and also
 its range. The naming convention of the parameters used in this software is
 slightly different than the naming convention used in [1]. Also, some of the
 default values of this software are different from the settings used in [1],
 so some care must be taken in setting these parameters.

 ``-mass (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)`` Specify the concentration
 parameter of the top level Dirichlet process in the BNPPhase model. This is
 the parameter ``alpha_0`` from [1]. The default value for ``mass`` is
 ``1e-3:1e2:1e3`` (i.e., initial value ``10``, minimum value ``1/1000`` and
 maximum value ``1000``).

 ``-varmass`` Specify the variance of the prior on the concentration parameter
 of the top level Dirichlet process. The default value for ``varmass`` is ``1``.
 
 ``-mu (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)`` Specify the concentration
 parameter of the bottom level Dirichlet processes in the BNPPhase model. This
 is the parameter ``alpha`` from [1]. The default value for ``mu`` is
 ``1e-3:1:1e3``.

 ``-varmu`` Specify the variance of the prior on the concentration parameter of
 the bottom level Dirichlet process. The default value for ``varmu`` is ``1``.

 ``-gamma (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)`` Specify the pseudo-
 observations for the beta/Bernoulli likelihood model. This is the parameter
 ``gamma`` from [1]. The default value for ``gamma`` is ``1e-3:1e-2:1e3``.

 ``-r (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)`` Specify the prior on the jump
 rate. This is the parameter ``r`` from [1]. The jump rate varies by location,
 but the initial jump rate will be set to the real number ``<INIT>`` at all
 locations for the initialisation of the MCMC.
 
 ``-b (<INIT>|<MIN>:<MAX>|<MIN>:<INIT>:<MAX>)`` Specify the range of the prior
 on the pseudo-observations for the minor allele frequency. This is the
 parameter ``beta`` from [1]. The default value for ``beta`` is
 ``1e-5:1e-3:1e-1``.

 ``-mass_prior <STRING>`` Specify the prior to use on the concentration of the
 top level Dirichlet process. Possible values for ``<STRING>`` are ``normal``,
 ``lognormal``, ``exponential`` or ``loguniform``. The default value is
 ``lognormal``. Note that if ``loguniform`` is used, then the ``varmass``
 parameter is ignored. The ``varmass`` parameter governs the variance of the
 prior before truncation to the range specified by the ``mass`` parameter

 ``-mu_prior <STRING>`` Specify the prior to use for the concentration of the
 bottom level Dirichlet process. Possible values for ``<STRING>`` are ``normal``,
 ``lognormal``, ``exponential`` or ``loguniform``. As for ``-mass_prior``,
 ``-varmu`` is ignored if ``<STRING>`` is ``loguniform``.

 ``-fixmass`` This flag specifies that the concentration parameter ``mass`` for
 the top level Dirichlet processes should be fixed, and not resampled during
 MCMC.  The default value for ``mass`` is ``1``. If a different value is desired,
 it can be provided with the ``-mass <INIT>`` parameter (and then the value of
 ``mass`` will be fixed to the real number ``<INIT>``). If ``-fixmass`` is
 provided, then the ``<MIN>:<MAX>`` or ``<MIN>:<INIT>:<MAX>`` forms for the flag
 ``-mass`` must not be used. The remaining flags (all beginning with ``-fix``)
 are of a similar form. The ``-fix`` flag can be prepended to the names of other
 parameters to institute a similar effect, as follows:

 ``-fixmu`` 

 ``-fixgamma`` 

 ``-fixr``

## REFERENCES:

[1] L.T. Elliott, Y.W. Teh. A nonparametric HMM for genetic imputation and
coalescent inference. 2016. Electronic Journal of Statistics.

[2] M. J. Daly, J. D. Rioux, S. F. Schaffner, T. J. Hudson, and R. S. Lander.
High-resolution haplotype structure in the human genome. 2001. Nature Genetics.


## CONTACT US:

If you would like to contact us for additional help, or to send us comments or
bug reports, please visit the ``https://github.com/BigBayes/BNPPhase`` website
and send us your message through github.


--------------------------------------------------------------------------------
