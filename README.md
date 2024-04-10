# [Facility location and pricing problem: Discretized mill price and exact algorithms](https://doi.org/10.1016/j.ejor.2022.11.052)


This archive is distributed under the [MIT License](LICENSE).

This repository contains supporting material for the paper 
 
    "Facility location and pricing problem: Discretized mill price and exact algorithms" by Yun Hui Lin, Qingyun Tian

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper.


## Cite

To cite the contents of this repository, please cite the paper

Below is the BibTex:

```
@article{lin2023facility,
  title={Facility location and pricing problem: Discretized mill price and exact algorithms},
  author={Lin, Yun Hui and Tian, Qingyun},
  journal={European Journal of Operational Research},
  volume={308},
  number={2},
  pages={568--580},
  year={2023},
  publisher={Elsevier}
}
```


## Description
The joint optimization of facility location and service charge arises in many industrial and business contexts. This paper investigates a facility location and mill pricing problem (FLMPr), where a company aims to maximize its profit by locating service facilities and setting appropriate service charges to customers. For each facility, the number of pricing levels are finite, and the company will select exactly one level for each facility if it is open. We visualize the problem from a fully decentralized perspective, i.e., each customer acts as an independent decision-maker. Under mill pricing, customers visiting the same facility encounter the same service charge. The problem is formulated as a bilevel program, in which the company makes location and pricing decisions at the upper level, and customers decide whether to seek the service from a facility at the lower level. To solve FLMPr, we leverage three types of closest assignment constraints to reformulate the problem as mixed-integer linear programs (MILPs), which can be directly solved by modern solvers. However, this approach suffers from a time-consuming solver compiling process and cannot handle large-scale instances effectively. This observation motivates us to design a branch-and-cut algorithm by exploring the bilevel structure and deriving a feasibility cut to efficiently eliminate bilevel infeasible solutions. Our extensive experiments reveal that the proposed algorithm can solve large-scale FLMPr satisfactorily and outperforms the MILP approach by a large margin. Finally, we conduct sensitivity analysis and draw interesting observations.

This repository provides data for the problem and code for the bilevel branch-and-cut algorithm.

## Replicating

- To run the code, you will need to make sure that you have already installed **Anaconda3**.

- You also need to install the solver **Gurobi** (license required)

## License

This software is released under the MIT license, which we report in file 'LICENSE'.
