# Graph-Drawing

This repository contains several graph drawing tools i made while working on my thesis for completing the BSc & MSc integrated diploma in Applied Mathematics and Physical Sciences @ National Technical University of Athens.  

      I will upload several parts of the project's code. For more details contact me at leandrosevag@hotmail.com

Monotone Drawings of Graphs is a new standard of graph drawing, introduced by Angelini et al. [[1]](#1).
A monotone drawing of a graph G is a straight-line drawing of G, s.t for every pair of nodes, there exists a path that is monotone in some direction (i.e. left-to-right). Several applications of monotone drawings includes optimal trajectory planning, map rotations and rearrangement, electronics design etc.



The main topics of interest is finding algorithms that reduce the size of the drawing. 

**Example  (from Angelini et al. paper) :**

![ex](https://i.ibb.co/x8fhRTT/angelini-et-al.png)

**Practical Example 1 (from our implementation):**

We can see how much better understanding of the graph we get from the monotone drawing. 

![mon](https://i.ibb.co/DzPNh1F/mon.png)

Monotone Drawings are based on the fact that humans have a specifiq process of analyzing visual information. Tt is shown that visualizations with high quality characteristics

**Practical Example 2 (from our implementation):**

Ternary tree with three levels, using three of the algorithms

![algos](https://i.ibb.co/qM00fps/algos.png)


## Technologies used in this project

![techs](https://i.ibb.co/qdfX2gc/techs.png)


## Introduction 

Every tree can be represented as a balanced set of open and closed parenthesis. The order of the nodes is also maintained. An example can be seen below. Starting from the root, every time we see an open parenthesis we move downwards and create a new node, else we go backwards. And we do this in counter-clockwise direction.

![trees](https://i.ibb.co/Xt9RTYZ/Untitled.png)

We know that every Tree can have a monotone drawing. In order to draw more general graphs, we compute a spanning-tree , draw it monotonically, and then add the extra edges. 

![mon](https://i.ibb.co/bLvL2ZF/Picture1.png)

This is why Monotone Drawings of Trees receeived much attention from researchers. 

In my thesis. I created a database from all distinct rooted trees with up to 15 nodes. I then implemented almost all of the known algorithms for monotone drawings of trees using Python and tested them on the database in order to analyze and compare their results in terms of grid size.
Also, i made a tool using HTML-JS-PHP-SQL-Python where a user can do one of the following tasks:

- Use and visualize monotone drawing algorithms for an input tree
- Create trees by entering their parethesis representation
- Import/Export graphml files
- Random Tree Generation
- Visual comparison of different algorithm grid-size utilization
- Access to the database (Retrieve and visualize data)

### The implemented algorithms with their grid size utilization can be seen below:

Algorithm papers:
  - Angelini et al. BFS and DFS [[1]](#1)
  - Kindermann et al. [[2]](#2)
  - He & He [[3]](#3)
  - Oikonomou & Symvonis  [[4]](#4)


![algorithms](https://i.ibb.co/gTRdn2r/algorithms.png)
      
      
      

Most of the algorithms use computational geometry and number theory for generating the "angles" assigned the edges of the graph. In particular, the Stern-Brocot Tree and the Farey Sequences are mostly used. The alogorithms of Oikonomou-Symvonis are the only ones not relying on number theory. They use a clever algorithmic way of assigning angles using basic geometry and they are the easiest to understand and implement. 

## Tree Generation

I include some algorithms for tree generation (as balanced set of parenthesis) and convertion to networkx datatype.

For a random tree we can compute a random valid string of parenthesis. The complexity of creating a random tree is   <img src="https://latex.codecogs.com/png.latex?\inline&space;O(n^{1.5})" title="O(n^{1.5})" />

For generating all distinct trees with n nodes, the computational complexity is   <img src="https://latex.codecogs.com/png.latex?\inline&space;O(\frac{4^{n}}{n^{1.5}})" title="O(\frac{4^{n}}{n^{1.5}})" />




## References

<a id="1">[1]</a> 
Journal of Graph Algorithms and Applications. Monotone Drawings of Graphs. 
Patrizio Angelini, Enrico Colasante, Giuseppe Di Battista,Fabrizio Frati,Maurizio Patrignani
[paper](https://www.emis.de/journals/JGAA/accepted/2012/Angelini+2012.16.1.pdf)

<a id="2">[2]</a> 
On Monotone Drawings of Trees
Philipp Kindermann, Andre Schulz, Joachim Spoerhase, Alexander Wolff
[paper](https://arxiv.org/pdf/1505.01410.pdf)

<a id="3">[3]</a> 
Optimal Monotone Drawings of Trees
Dayu He, Xin He
[paper](https://arxiv.org/pdf/1604.03921v1.pdf)

<a id="4">[4]</a> 
SIMPLE COMPACT MONOTONE TREE DRAWINGS 
ANARGYROS OIKONOMOU, ANTONIOS SYMVONIS
[paper](https://arxiv.org/pdf/1708.09653.pdf)







