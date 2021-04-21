# MOLSD
Multi-objective Iterated Local Search Algorithm based on Decomposition (MOLSD)

The algorithm called MOLSD has been designed to solve a wide range of Multi-objective Permutation Flow-shop Scheduling Problems (MoPFSP). The MoPFSP represents a set of real-world scheduling problems from the industry and engineering. The problem has a scalable number of jobs that have to be processed by a set of machines sequentially. The goal is to optimize one or multiple objective functions, commonly related to the completion time. 

MOLSD has been designed in C++ using the concepts of Object-oriented programming. We have used the Eclipse environment on Ubuntu. The algorithm is based on the traditional MOEA/D framework, which decomposes a multi-objective problem into a set of N scalar subproblems. Each subproblem is associated with a weighted vector. All the subproblems are optimized simultaneously in a collaborative manner using a concept of neighborhood between them. Each subproblem applies a search-based operator (such as crossover and mutation from the genetic operators, probabilistic models, and local search). 

In the paper "A Decomposition-based Local Search Algorithm for Multi-objective Sequence Dependent Setup Times Permutation Flowshop Scheduling" published on "2018 IEEE Congress on Evolutionary Computation", we present 1) the designed MOLSD, the 2) problem targeted, and 3) the experimental studies conducted.

The files present in the repository show the source code in which includes several classes. The algorithm has been designed to have a general framework, robust and flexible.

The published paper is available at https://ieeexplore.ieee.org/document/8477711 or you can download it here https://www.researchgate.net/publication/325020567_A_Decomposition-Based_Local_Search_Algorithm_for_Multi-Objective_Sequence_Dependent_Setup_Times_Permutation_Flowshop_Scheduling

