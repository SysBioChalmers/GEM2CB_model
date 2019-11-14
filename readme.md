# GEM2CB
Some functions to bulid cybernetic model from GEMs
Updated: 2019-10-27

Requirements,

- [cobrapy](https://opencobra.github.io/cobrapy/)
- [cvxpy](https://github.com/cvxgrp/cvxpy)
- [IBM Cplex](https://www.ibm.com/se-en/analytics/cplex-optimizer)

# Notice <br />
Please add `./ComplementaryScripts/My_def` to sys.path to ues functions in My_def   <br />
    - `import sys`    <br />
    - `sys.path.extend(['./ComplementaryScripts/My_def'])` 
    

# Structre,
Published pipeline:

Reduced GEMS --> EFMs --> MYA --> Parameter fitting --> Cybernetic model

Current pipeline:

**GEMs --> Yield space --> MYA** --> Parameter fitting --> Cybernetic model

Important Scrips/ Functions introductions:

GEM2pathways.py: 
- get yield space from GEMs,
- core function : opt yield (yield = p/s as opjective)
- need substract and productions list.

ConvexHull_yield.py:
- MYA 
- Hull all --> Hull cutoff
- Experiment data: (Bugs?? test again )


 
Parameter fitting & ODE process 
- TODO 


# Contact Us

lhao@chalmers.se



