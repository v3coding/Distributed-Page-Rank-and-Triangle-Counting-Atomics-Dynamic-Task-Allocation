Modified the Page Rank and Triangle Counting Algorithms to work on a distributed system using Barriers for synchronization


For Page Rank :

in : page_rank_parallel.cpp

All of the parallel implementations are written entirely by me, modified from and built upon the Serial implementation

The Page Rank Algorithm calculates the rank of a Node based on how many other Nodes "point" to it, and their respective ranks.

Three types of algorithmic methods were tested -- 

Edge Based Decomposition : 

    Task decomposition is done by first breaking the task down in to the number of edges on the Graph 

    The number of Edges, E, is divided by the number of Processes, P

    The number of edges to be computed by each Process is then divided approximately evenly, with the last process taking up the remainder of the division

    Then, we count the number of edges per Node until we get start and end nodes for each process on the Graph

    We then pass the start and edge nodes to the page rank parallel method and it computes the page rank in parallel 

Dynamic Decomposition : 

    A main thread keeps track of and serves the next Node to be processed on graph G to the next available process

    Each process works in paralell until the last Node is processesed 

Vertex Based Decomposition : 

    The Graph is split up in to approximately V/P vertexes per process with the last process taking up the remainder

    The processes work in parallel to compute the Page Rank 



For Triangle Counting : 

in : triangle_counting_parallel.cpp

All of the parallel implementations are written entirely by me, modified from and built upon the Serial Implementation

The Triangle Counting Algorithm counts how many Nodes on a Graph form a triangle with Edges connecting them

Three types of algorithmic methods were tested -- 

Edge Based Decomposition : 

    Analogous to the case above for Page Rank

Dynamic Decomposition 

    Analogous to the case above for Page Rank

Vertex Based Decomposition :

    Analogous to the case above for Page Rank

    





