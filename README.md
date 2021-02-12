# Private-Stream-Aggregation-algorithm
Realised end 2020 for a semester project at EURECOM.

## Presentation
This project aims at implementing the Private Stream Aggregation algorithm proposed by Elaine Shi, T-H. Hubert Chan, Eleanor Rieffel, Richard Chow and Dawn Song (www.elaineshi.com/docs/ndss2011.pdf) using *Charm Crypto*.
The goal of this scheme is basically to allow N users to engage in statistical operations over their shared data without revealing their individual data. An aggregator, which could be a server, is able to do the sum of all users'data without compromising each individual privacy. Hence, each users and the aggregator use different keys (contrary to homomorphic encryption) and this means that users neither need to trust each other nor the aggregator.
The scheme provides also distributed differential privacy to protect against indirect violation of privacy (e.g compromised participants ally with the aggregator to infer data from other participants).

## Content
* PSA_charm_scheme.py: main file with the implementation of the algorithm and the definition of parameters.
* Pollard_algorithm.py: implementation of Pollard algorithm for discrete logarithm computation adapted for Charm Crypto.
* DD_privacy: generates graph to compare the distributed differential privacy scheme against a naive scheme.

It is really recommended to read the paper to understand the implementation.

## Dependencies
* Charm Crypto

Charm Crypto is a Python framework made for prototyping complex cryptosystems, promoting the reuse of components. You can follow this link to install it https://jhuisi.github.io/charm/install_source.html.
We mainly use it for generating and using Shnorr group.
* Scipy

Scipy ecosystem gives you a lot of tools for computation, data management, mathematical experimentation. It is included in most Python distributions, follow this link for installation information https://www.scipy.org/install.html. We only import scipy.stat to build the two-sided geometric distribution.

## Use
Run the PSA_charm_scheme after defining the desired parameters:
* bits : number of bits used for prime numbers
* delta : limit value of user's data (user's value included in {0,1,...,delta})
* gamma : proportion of uncompromised paraticipants
* nb_users : number of users
* eps, d : (eps, d)-distributed-differential-privacy parameters with delta >= eps/3 and 0 < d < 1

(Do not hesitate to look at the code since it is fully commented ! If you have a question about the code you can send an email to lardy@eurecom.fr)





