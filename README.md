# CorrLDA

The C implementation for Correspondence latent Dirichlet allocation with collapsed gibbs sampling estimation.

[Blei and Jordan. 2003] http://www.cs.columbia.edu/~blei/papers/BleiJordan2003.pdf

Fixed-point iteration method are used for hyper-parameter update. [Minka 2000] http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf


## Build

* Requirements : The GNU Scientific Library (GSL)

* Edit the Makefile appropriately and type "make" command.


## Usage

asymcorrlda [-I maxiter] [-K n_classes] [-A alpha] [-B beta] [-G gamma] [-S random_seed] doc tag model


-I (int) : The number of Gibbs iterations.

-K (int) : The number of topics.

-A (float) : The initial settings of an asymmetric Dirichlet prior over the document-topic distributions. These values are optimized during the Gibbs iterations.

-B (float) : The initial setting of an symmetric Dirichlet prior over the topic-word distributions. This value are optimized during the Gibbs iterations.

-G (float) : The initial setting of an symmetric Dirichlet prior over the topic-tag distributions. This value are optimized during the Gibbs iterations.

-S (int) : a seed value for the random number generator (to fix the initial setting of the Gibbs iteration).

doc : An input file for word-counts of documents.

tag : An input file for tag-counts of documents.

model : a prefix for output model files.


## A format for input files (doc and tag)

Each line indicates each document.

For each document, counts of each unique word are written as "(WordID):(Count)", where WordID is an index for words in dataset (1-origin, so the first index of words are 1, not zero), and Count is the number of times the word appeared on the document. Multiple values are separated by a space.


E.g.)

12:1 353:1 416:3 636:2 670:1 713:1

38:2 72:2 109:2 265:1

18:1 40:1 98:1 251:1 265:3 411:1 743:1

...


## Output files

* [model].theta : Topic distributions for each document. (The number of documents) x (The number of topics) matrix.

* [model].phi : Word distributions for each topic. (The size of word vocabulary) x (The number of topics) matrix.

* [model].psi : Tag distributions for each topic. (The size of tag vocabulary) x (The number of topics) matrix.

* [model].hyper : Optimization processes for hyper parameters. Each line shows parameters at each step of Gibbs iterations. From left, alpha parameters (for each topic), a beta parameter, and a gamma parameter.

* [model].lik : Log-likelihood of the model in each step of iterations.
