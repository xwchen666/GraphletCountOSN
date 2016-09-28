# Mining Graphlet Counts in Online Social Networks

We propose an efficient random walk-based framework to estimate the subgraph
counts.  Our framework generates samples by leveraging consecutive steps of the
random walk as well as by observing neighbors of visited nodes. Using the
importance sampling technique, we derive unbiased estimators of the subgraph
counts. To make better use of the degree information of visited nodes, we also
design an improved estimator, which increases the efficiency of the estimate at
no additional cost. We conduct extensive experi- mental evaluation on
real-world OSNs to confirm our theoretical claims. The experiment results show
that our estimators are unbiased, accurate, efficient and better than the
state-of-the-art algorithm. For the Weibo graph with more than 58 million
nodes, our method produces estimate of triangle count with an error less than
0.05 using only 20 thousands sampled nodes. Detailed comparison with the
state-of-the-art method demonstrates that our algorithm is 4 to 5 times more
accurate.


## How to run

- mkdir build
- cd build
- cmake ..
- make


## Datasets

We list the detailed information of the LCCs of the datasets used in our paper.

| **Name**     | **Node**  | **Edges** | **Description** | **Link**                                                          |
|----------|-------|-------|----------------------------------------------------------------|---------------------------------------------------------------|
| Epinion  | 76K   | 406K  | Trust network from the onlinesocial network Epinion            | http://snap.stanford.edu/data/soc-Epinions1.html |
| Slashdot | 77K   | 469K  | Friend/foe links between the usersof Slashdot social network   | http://snap.stanford.edu/data/soc-Slashdot0811.html |
| Pokec    | 1.6M  | 22.3M | Friendship network from theSlovak social network Pokec         | http://konect.uni-koblenz.de/networks/soc-pokec-relationships |
| Flickr   | 2.2M  | 22.7M | Social network of Flickr users andtheir friendship connections | http://konect.uni-koblenz.de/networks/flickr-growth           |
| Twitter  | 21.3M | 265M  | Graph about who follows whomon Twitter | http://networkrepository.com/soc-twitter-2010.php             |
| Weibo    | 58.7M | 261M  | A micro-blogging service withmillions of users in China.       | http://networkrepository.com/soc-sinaweibo.php                |

## Usage

Our implementations of the proposed algorithm is in *IMPR.h* and *IMPR.cpp*. Please refer to our paper for detailed derivation of the unbiased estimators. 

## Details

We also include the algorithms from the following three papers

- Ahmed, Nesreen K., et al. "**Efficient graphlet counting for large networks.**", ICDM'15
- Hočevar, Tomaž, and Janez Demšar. "**A combinatorial approach to graphlet counting.**" Bioinformatics 30.4 (2014): 559-565.  
- Wang, Pinghui, et al. "**Efficiently estimating motif statistics of large networks.**" TKDD'14 


## References

- Xiaowei Chen, John C.S. Lui, **Mining Graphlet Counts in Online Social Networks**. In ICDM'16, to appear

